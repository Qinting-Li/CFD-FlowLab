#include "DSMCSolver.h"
#include "DSMCMolecule.h"
#include "DSMCConsts.h"
#include "DSMCGrid.h"
#include "DSMCSample.h"
#include "DSMCIndex.h"
#include "DSMCCollision.h"
#include "TK_Parse.h"
#include "TK_Time.h"
#include "Region.h"
#include <cmath>
using namespace std;

namespace PHSPACE
{

DSMC::DSMC()
{
    region = new Region();
    testTime = new TimeSpan();
    species = NULL;
    SPM = NULL;

    dsmcgrid = NULL;
    molecule = NULL;
    sample = NULL;
    index = NULL;
    collision = NULL;
}

DSMC::~DSMC()
{
    delete region; region = NULL;
    delete testTime; testTime = NULL;
    delete [] species;  species = NULL;
    DelPointer3 (SPM) ; SPM = NULL;

    int nZones = PHMPI::GetNumberofGlobalZones();
    for(int izone = 0; izone < nZones; ++izone)
    {
        delete dsmcgrid[izone];
        delete molecule[izone];
        delete sample[izone];
        delete index[izone];
        delete collision[izone];
    }
    delete [] dsmcgrid; dsmcgrid = NULL;
    delete [] molecule; molecule = NULL;
    delete [] sample; sample = NULL;
    delete [] index; index = NULL;
    delete [] collision; collision = NULL;
}

void DSMC::Run()
{
    using namespace PHMPI;
    double last_time[3];
    TIME_SPACE::StartTime(last_time);
    cout<<"Begin DSMC simulation"<<endl;

    PreProcess();
    PrintToWindow("PreProcess completed \n");
    int maxSimuStep = GlobalDataBase::GetIntParaFromDB("maxSimuStep");
    int nsample = GlobalDataBase::GetIntParaFromDB("nsample");
    int nout = GlobalDataBase::GetIntParaFromDB("nout");
    int intervalStepFlow = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");

    for (int nstep = 0; nstep < maxSimuStep; ++ nstep)
    {
        PrintToWindow("nstep = ",nstep,"  ");
        cout<<dsmcgrid[0]->GetnMolecule()<<endl;
        if(nstep == nsample) 
            InitSample();

        MoleculeMove();

        MoleculeIndex();

        MoleculeCollision();

        FieldSample();

        if (nstep % nout == 0 || nstep == maxSimuStep)
            PostProcess();

        if (nstep % intervalStepFlow == 0 || nstep == maxSimuStep)
            WriteResumeFiles();
    }

    //PrintToWindow("DeleteMemory \n");
    //DeleteMemory();
    //数字转字符串的方法，lijin，2019-2-25
    //begin
    //int num=1;
    //RDouble a=123.33;
    //stringstream filename;
    //string filename2;
    //filename<<a;
    //filename>>filename2;
    //cout<<"filename"<<filename2<<endl;
    //end

}
void DSMC::PostProcess()
{
    RDouble SUMM[13];
    int nZones = PHMPI::GetNumberofGlobalZones();
    for(int izone = 0; izone < nZones; ++izone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(izone);

        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));

        int nTotalCell = grid->GetNTotalCell();
        RDouble * cellVolume = grid->GetCellVolume();
        RDouble *** CS = sample[izone]->GetCS();
        RDouble ** VAR = sample[izone]->GetVAR();
        RDouble *** VARSP = sample[izone]->GetVARSP();
        RDouble ** CSR = sample[izone]->GetCSR();

        int K;
        for(int icell = 0; icell < nTotalCell; ++ icell)
        {
            RDouble A = FNUMI[izone]/cellVolume[icell]/NSMP;
            for(K = 0; K < 13; ++K)
                SUMM[K] = 0.0;
            for(K = 0; K < 21; ++K)
                VAR[K][icell] = 0.0;

            for(int isp = 0; isp < nSpecies; ++isp)
            {
                SUMM[0] += CS[0][icell][isp];
                SUMM[1] += CS[1][icell][isp];
                SUMM[2] += CS[1][icell][isp]*species[isp].mass;
                for(K = 1; K <= 3 ; ++K)
                {
                    SUMM[K+2] += CS[K+1][icell][isp]*species[isp].mass;
                    if(CS[1][icell][isp] > 0.1)
                    {
                        VARSP[K+1][icell][isp] = CS[K+4][icell][isp]/CS[1][icell][isp];
                        VARSP[K+8][icell][isp] = CS[K+1][icell][isp]/CS[1][icell][isp];
                    }
                }

                SUMM[6] += (CS[5][icell][isp] + CS[6][icell][isp] + CS[7][icell][isp])*species[isp].mass;

                if(CS[1][icell][isp] > 0.5)
                    SUMM[7] += CS[5][icell][isp] + CS[6][icell][isp] + CS[7][icell][isp];

                if(species[isp].rotdof > 0)
                {
                    SUMM[8] += CSR[icell][isp];
                    SUMM[9] += CS[1][icell][isp]*species[isp].rotdof;
                }

                SUMM[10] += CS[5][icell][isp]*species[isp].mass;
                SUMM[11] += CS[6][icell][isp]*species[isp].mass;
                SUMM[12] += CS[7][icell][isp]*species[isp].mass;
            }

            for(int isp = 0; isp < nSpecies; ++isp)
            {
                VARSP[0][icell][isp] = CS[1][icell][isp];
                VARSP[1][icell][isp] = 0.0;
                VARSP[6][icell][isp] = 0.0;
                VARSP[7][icell][isp] = 0.0;
                VARSP[8][icell][isp] = 0.0;

                if(SUMM[1] > 0.1)
                {
                    VARSP[1][icell][isp] = CS[1][icell][isp]*A;
                    if(species[isp].rotdof > 0 &&  CS[1][icell][isp] > 0.5)
                        VARSP[6][icell][isp] = (2.0/BOLTZ)*CSR[icell][isp]/species[isp].rotdof/CS[1][icell][isp];
                }
                VARSP[5][icell][isp] = 0.0;
                for(K = 1; K <= 3 ; ++K)
                {
                    VARSP[K+1][icell][isp] = (species[isp].mass/BOLTZ)*(VARSP[K+1][icell][isp]-VARSP[K+8][icell][isp]*VARSP[K+8][icell][isp]);
                    VARSP[5][icell][isp] += VARSP[K+1][icell][isp];
                }
                VARSP[5][icell][isp] /= 3.0;
            }

            VAR[1][icell] = SUMM[1];
            VAR[2][icell] = SUMM[0];

            if(SUMM[1] > 0.5)
            {
                VAR[3][icell] = SUMM[1]*A;
                VAR[4][icell] = VAR[3][icell]*SUMM[2]/SUMM[1];
                VAR[5][icell] = SUMM[3]/SUMM[2];
                VAR[6][icell] = SUMM[4]/SUMM[2];
                VAR[7][icell] = SUMM[5]/SUMM[2];
                RDouble UU = VAR[5][icell]*VAR[5][icell] +  VAR[6][icell]*VAR[6][icell] + VAR[7][icell]*VAR[7][icell];
                if(SUMM[1] >1) 
                {
                    VAR[8][icell] = abs(SUMM[6] - SUMM[2]*UU)/(3.0*BOLTZ*SUMM[1]);
                    VAR[19][icell] = abs(SUMM[10] - SUMM[2]*VAR[5][icell]*VAR[5][icell])/(BOLTZ*SUMM[1]);
                    VAR[20][icell] = abs(SUMM[11] - SUMM[2]*VAR[6][icell]*VAR[6][icell])/(BOLTZ*SUMM[1]);
                    VAR[21][icell] = abs(SUMM[12] - SUMM[2]*VAR[7][icell]*VAR[7][icell])/(BOLTZ*SUMM[1]);
                }
                else
                {
                    VAR[8][icell] = 1.0;
                    VAR[19][icell] = 1.0;
                    VAR[20][icell] = 1.0;
                    VAR[21][icell] = 1.0;
                }

                if(SUMM[9] > 0.01)
                    VAR[9][icell] = (2.0/BOLTZ)*SUMM[8]/SUMM[9];
                else
                    VAR[9][icell] = 0.0;

                VAR[10][icell] = 0.0;

                RDouble DOF = 3.0 + SUMM[9]/SUMM[1];

                VAR[11][icell] = (3*VAR[8][icell]+SUMM[9]/SUMM[1]*VAR[9][icell])/DOF;
                VAR[18][icell] = VAR[3][icell]*BOLTZ*VAR[8][icell]; 

            }


        
        }
    
    }


}

void DSMC::FieldSample()
{
    NSMP += 1;
    int nZones = PHMPI::GetNumberofGlobalZones();
    int nweight = Getnweight();
    int MMVM = GetMMVM();
    int MMRM = GlobalDataBase::GetIntParaFromDB("MMRM");
    RDouble wfactor;
    for(int izone = 0; izone < nZones; ++izone)
    {
        int nMolecule = dsmcgrid[izone]->GetnMolecule();
        int * moleculeLocation = molecule[izone]->GetMoleculeLocation();
        for(int imolecule = 0; imolecule < nMolecule; ++ imolecule)
        {
            int icell = moleculeLocation[imolecule];
            if(nweight == 0)
            {
                wfactor = 1.0;
            }
            else
            {
                wfactor = dsmcgrid[izone]->GetCellWeight()[icell];
            }

            int isp = molecule[izone]->GetMoleculeSpecies()[imolecule];
            sample[izone]->GetCS()[0][icell][isp] += 1.0;
            sample[izone]->GetCS()[1][icell][isp] += wfactor;

            for(int k = 0; k < 3; ++ k)
            {
                sample[izone]->GetCS()[k+2][icell][isp] += wfactor*molecule[izone]->GetMoleculeVelocity()[k][imolecule]; 
                sample[izone]->GetCS()[k+5][icell][isp] += wfactor*molecule[izone]->GetMoleculeVelocity()[k][imolecule]*molecule[izone]->GetMoleculeVelocity()[k][imolecule]; 
            }

            if(MMRM > 0)
                sample[izone]->GetCSR()[icell][isp] += wfactor*molecule[izone]->GetMoleculeRotationalEnergy()[imolecule];

            if(MMVM > 0)
                sample[izone]->GetCSV()[0][icell][isp] += wfactor*molecule[izone]->GetMoleculeVibrationLevel()[imolecule]*BOLTZ*species[isp].vibtemp;
        
        }
    }


}
void DSMC::PreProcess()
{   
    using namespace PHMPI;

    ReadSpecies();    

    //ZoneToCPU();

    ReadGeometry();    

    UpdateAllZonesPara();

    InitGeometry();  

    PrintToWindow("AllocateMemory \n");
    AllocateMemory();

    PrintToWindow("Calculating reference cell volume \n");
    CalReferenceCellVolume();

    int nZones = GetNumberofGlobalZones();
    for(int izone = 0; izone < nZones; ++izone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(izone);
        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));

        //PrintToWindow("Calculating cell weight for izone=",izone," \n");
        dsmcgrid[izone]->CalCellWeight(grid,GetReferenceCellVolume());

    }

    PrintToWindow("Calculating time step \n");
    CalTimeStep();

    PrintToWindow("Calculating the ratio of real molecule number represented by a simulator \n");
    CalRealtoSimRatio();

    int nstart = GlobalDataBase::GetIntParaFromDB("nstart");

    if(nstart == 0)
    {
        InitData();

        MoleculeIndex();

        InitSample();

    }
    else if(nstart == 1)
    {
        
        ReadRestart();

    }
    else 
    {
        TK_Exit::UnexpectedVarValue("nstart = ",  nstart);
    }

    

}

void DSMC::ReadGeometry()
{
    using namespace PHMPI;

    //TimeSpan * testTime = new TimeSpan();

    region->ReadGrid();

    region->UpdateAllZonesPara();

    testTime->ShowTimeSpanToLogFile("Grid Reading");

    int nZones = GetNumberofGlobalZones();
    dsmcgrid = new DSMCUnstructGrid * [nZones];

    for (int izone = 0; izone < nZones; ++ izone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(izone);

        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));
        dsmcgrid[izone] = new DSMCUnstructGrid();

        int nTotalNode = grid->GetNTotalNode();

        RDouble *x = grid->GetX();
        RDouble *y = grid->GetY();
        RDouble *z = grid->GetZ();

        RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

        for(int inode = 0; inode < nTotalNode; ++ inode)
        {
            x[inode] *= gridScaleFactor;
            y[inode] *= gridScaleFactor;
            z[inode] *= gridScaleFactor;
        }

        grid->SetX(x);
        grid->SetY(y);
        grid->SetZ(z);

        //int **cell2face  = grid->GetCell2Face();
        //int * face_number_of_each_cell = grid->GetFaceNumberOfEachCell();

        //for (int iface = 0; iface < face_number_of_each_cell[0]; ++ iface)
        //{
        //  int face = cell2face[0][iface];
        //}

        //int * cell2node = grid->GetCell2Node(); 

        //RDouble * xfc = grid->GetFaceCenterX();
        //RDouble * yfc = grid->GetFaceCenterY();
        //RDouble * zfc = grid->GetFaceCenterZ();

        //RDouble * xfn = grid->GetFaceNormalX();
        //RDouble * yfn = grid->GetFaceNormalY();
        //RDouble * zfn = grid->GetFaceNormalZ();

        //RDouble * xcc = grid->GetCellCenterX();
        //RDouble * ycc = grid->GetCellCenterY();
        //RDouble * zcc = grid->GetCellCenterZ();
        //
        //
        //RDouble * area = grid->GetFaceArea();

        //int * node_number_of_each_face = grid->GetNodeNumberOfEachFace();

        //int * face2node = grid->GetFace2Node();
        //int * left_cell_of_face = grid->GetLeftCellOfFace();
        //int * right_cell_of_face = grid->GetRightCellOfFace();

        //CompositeUnstructBC **bcr = grid->GetBCRecord();
        //ostringstream oss;
        //for(int iface = 0; iface < nTotalFace; ++ iface)
        //{
        //  int le = left_cell_of_face[iface];
        //  int re = right_cell_of_face[iface];
        //  oss<<le<<" "<<re<<endl;
        //}
        //oss<<"========================"<<endl;
        //for(int iface = 0; iface < nBoundFace; ++ iface)
        //{
        //  int le = left_cell_of_face[iface];
        //  int re = right_cell_of_face[iface];
        //  oss<<le<<" "<<re<<" "<<bcr[iface]->GetKey()<<endl;
        //
        //}
        //WriteLogFile(oss.str());
        //cin.get();

        dsmcgrid[izone]->BoundaryCount(grid);

    }

}



void DSMC::UpdateAllZonesPara()
{
    this->region->UpdateAllZonesPara();
}

void DSMC::InitGeometry()
{
    this->region->InitGeometry();
}

void DSMC::AllocateMemory()
{
    int nZones = PHMPI::GetNumberofGlobalZones();

    molecule = new Molecule * [nZones];
    sample = new Sample * [nZones];
    index = new Index *[nZones];
    collision = new Collision *[nZones];

    for(int izone = 0; izone < nZones; ++izone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(izone);
        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));
        int nTotalCell = grid->GetNTotalCell();

        molecule[izone] = new Molecule();
        sample[izone] = new Sample();
        index[izone] = new Index();
        collision[izone] = new Collision(nTotalCell);
    }

    for (int izone = 0; izone < nZones; ++ izone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(izone);

        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));

        int nTotalCell = grid->GetNTotalCell();

        index[izone]->SetStartAddressOfMoleculeInCells(nTotalCell);
        index[izone]->SetNumberOfMoleculeInCells(nTotalCell);
        index[izone]->SetMoleculeNoInOrderOfCells(MaxMoleculeNumber);

        sample[izone]->SetCSR(nTotalCell,nSpecies);

    }
}

void DSMC::DeleteMemory()
{
    int nZones = PHMPI::GetNumberofGlobalZones();

    for (int izone = 0; izone < nZones; ++ izone)
    {
        int * startAddressOfMoleculeInCells = index[izone]->GetStartAddressOfMoleculeInCells();
        int * numberOfMoleculeInCells = index[izone]->GetNumberOfMoleculeInCells();
        int * moleculeNoInOrderOfCells = index[izone]->GetMoleculeNoInOrderOfCells();

        delete [] startAddressOfMoleculeInCells ;
        delete [] numberOfMoleculeInCells;
        delete [] moleculeNoInOrderOfCells;
    }
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 函数头部说明：                                                             
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:   例如：CalRealtoSimRatio()                      
// AUTHOR(S):  Li Jin                                                          
// PURPOSE:    计算一个模拟分子所代表的真实分子                                                  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DSMC::CalRealtoSimRatio()
{
    using namespace PHMPI;

    int nSpecies = GetnSpecies();
    Species * species = this->species;

    int nZones = GetNumberofGlobalZones();
    FNUMI.resize (nZones);
    DTMI.resize (nZones);
    TIME.resize (nZones);
    TIMI.resize (nZones);
    Zone2CPU.resize(nZones);

    int nsim = GlobalDataBase::GetIntParaFromDB("nsim");
    RDouble inflowNumberDensity = GlobalDataBase::GetDoubleParaFromDB("inflowNumberDensity"); 

    ostringstream oss;
    oss << setiosflags(ios::right);
    oss << setprecision(3);
    //oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);
    int wordwidth = 4;
    oss << "Simulators for each specie:" << endl;

    for(int i = 0; i < nSpecies; ++i)
    {
        oss << setw(wordwidth) << species[i].id <<" " << nsim*species[i].percent << "\n";
    }
    WriteLogFile(oss.str());

    RDouble weight = 1.0;

    for(int izone = 0; izone < nZones; ++izone)
    {
        FNUMI[izone] = inflowNumberDensity*GetReferenceCellVolume()/weight/nsim;
    }

    for(int izone = 0; izone < nZones; ++izone)
    {
        DTMI[izone] = timeStep;
    }

}



void DSMC::CPUtoZone()
{
    using namespace PHMPI;
    
    int nZones = GetNumberofGlobalZones();
    int nCPU = GetNumberOfProcessor();

    int nS,nS0,Mflag;
    Mflag = 0;
    for(int i = 0; i < nZones ; ++i)
    {
        int j = ZoneNo[i];
        for(int n = 0; n < nCPU; ++n)
        {
            if(n < ZoneLeft)
            {
                nS = ZonePerCPU+1;
                nS0 = n*nS;
            }
            else
            {
                nS = ZonePerCPU;
                nS0= ZoneLeft*(nS+1)+(n-ZoneLeft)*nS;
            }
            if(i > nS0 && i <= nS0+nS) 
            {
                Mflag = j;
                exit(0);
            }
        
        }
        //if(Mflag < 0) call Exit(-5201);
        Zone2CPU[j] = Mflag;
    }



}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 函数头部说明：                                                              +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:   例如：CalTimeStep(Molecule *& molecule)                         +
// AUTHOR(S):  Li Jin                                                          +
// PURPOSE:    计算时间步长                                                    +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DSMC::CalTimeStep()
{
    int nSpecies = GetnSpecies();
    Species * species = this->species;
    
    RDouble mixtureMass = 0.0;
    for(int isp = 0; isp < nSpecies; isp++)
    {
        mixtureMass += species[isp].mass * species[isp].percent;
    }

    RDouble inflowDensity = GlobalDataBase::GetDoubleParaFromDB("inflowDensity");
    RDouble inflowTemperature = GlobalDataBase::GetDoubleParaFromDB("inflowTemperature");
    RDouble inflowVelocity = GlobalDataBase::GetDoubleParaFromDB("inflowVelocity");
    RDouble inflowMachNumber = GlobalDataBase::GetDoubleParaFromDB("inflowMachNumber");

    RDouble inflowNumberDensity = inflowDensity/mixtureMass;
    GlobalDataBase::UpdateData("inflowNumberDensity", &inflowNumberDensity, PHDOUBLE, 1);

    RDouble specificHeatRatioAtConstantVolume=0.0;
    RDouble specificHeatRatioAtConstantPressure=0.0;
    for(int isp = 0; isp < nSpecies; isp++)
    {
        specificHeatRatioAtConstantVolume += 0.5*(species[isp].rotdof + 3)*BOLTZ*inflowNumberDensity*species[isp].percent/inflowDensity;
        specificHeatRatioAtConstantPressure += 0.5*(species[isp].rotdof + 5)*BOLTZ*inflowNumberDensity*species[isp].percent/inflowDensity;
    }

    RDouble specificHeatRatio = specificHeatRatioAtConstantPressure/specificHeatRatioAtConstantVolume;
    RDouble gasConstant = BOLTZ/mixtureMass;

    GlobalDataBase::UpdateData("specificHeatRatio", &specificHeatRatio, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("gasConstant", &gasConstant, PHDOUBLE, 1);
    
    RDouble inflowSonicSpeed = sqrt(specificHeatRatio*gasConstant*inflowTemperature);
    GlobalDataBase::UpdateData("inflowSonicSpeed", &inflowSonicSpeed, PHDOUBLE, 1);

    if(inflowVelocity < 0 && inflowMachNumber > 0) 
    {
        inflowVelocity = inflowSonicSpeed*inflowMachNumber;
    }
    else if(inflowVelocity > 0 && inflowMachNumber < 0)
    {
        inflowMachNumber = inflowVelocity/inflowSonicSpeed;
    }

    ostringstream oss;
    oss<<"specificHeatRatio："<<specificHeatRatio<<endl;
    oss<<"inflowSonicSpeed："<<inflowSonicSpeed<<endl;
    oss<<"inflowVelocity："<<inflowVelocity<<endl;
    oss<<"gasConstant："<<gasConstant<<endl;
    oss<<"inflowNumberDensity："<<inflowNumberDensity<<endl;

    SPM = NewPointer3<RDouble>(7, nSpecies,nSpecies);

    for(int isp = 0; isp < nSpecies; isp++)
    {
        for(int jsp = 0; jsp < nSpecies; jsp++)
        {
            SPM[0][isp][jsp] = 0.25*PI*pow(species[isp].diam+species[jsp].diam,2);
            SPM[1][isp][jsp] = 0.5*(species[isp].tref+species[jsp].tref);
            SPM[2][isp][jsp] = 0.5*(species[isp].omega+species[jsp].omega);
            SPM[3][isp][jsp] = 0.5*(species[isp].alpha+species[jsp].alpha);
            SPM[4][isp][jsp] = species[isp].mass*species[jsp].mass/(species[isp].mass+species[jsp].mass);
            SPM[5][isp][jsp] = gamma(static_cast <RDouble> (2.5-SPM[2][isp][jsp]));
            if(species[isp].rotdof > 0 && species[jsp].rotdof > 0)
            {
                SPM[6][isp][jsp] = 0.5*(species[isp].rotc1 + species[jsp].rotc1);//转动碰撞数，目前是读入的是倒数，一般是5，现在读入的是0.2
            }
            else if(species[isp].rotdof > 0 && species[jsp].rotdof == 0)
            {
                SPM[6][isp][jsp] = species[isp].rotc1;
            }
            else if(species[isp].rotdof == 0 && species[jsp].rotdof > 0)
            {
                SPM[6][isp][jsp] = species[jsp].rotc1;
            }
        }
    }

    RDouble inflowMeanFreePath = 0.0;
    RDouble inflowMeanCollisionTime = 0.0;
    RDouble inflowMeanCollisionRate = 0.0;

    for(int isp = 0; isp < nSpecies; isp++)
    {
        RDouble ipercent = species[isp].percent;
        RDouble MeanCollisionRate = 0.0;
        RDouble MeanFreePath = 0.0;

        for(int jsp = 0; jsp < nSpecies; jsp++)
        {
            RDouble dref = sqrt(SPM[0][isp][jsp]/PI);
            RDouble tref = SPM[1][isp][jsp];
            RDouble omega = SPM[2][isp][jsp];
            RDouble Mref = SPM[4][isp][jsp];
            RDouble jpercent = species[jsp].percent;
            RDouble nq = inflowNumberDensity*jpercent;
            MeanCollisionRate += 2*sqrt(PI)*dref*dref*nq*pow(inflowTemperature/tref,1-omega)*sqrt(2*BOLTZ*tref/Mref);
            MeanFreePath += PI*dref*dref*nq*pow(tref/inflowTemperature,omega-0.5)*sqrt(1+species[isp].mass/species[jsp].mass);
        }
        inflowMeanCollisionRate += MeanCollisionRate*ipercent;
        inflowMeanFreePath += (1.0/MeanFreePath)*ipercent;
    }
    inflowMeanCollisionTime = 1.0/inflowMeanCollisionRate;

    GlobalDataBase::UpdateData("inflowMeanFreePath", &inflowMeanFreePath, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("inflowMeanCollisionTime", &inflowMeanFreePath, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("inflowMeanCollisionRate", &inflowMeanFreePath, PHDOUBLE, 1);

    RDouble forceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble inflowKnudsenNumber = inflowMeanFreePath/forceReferenceLength;
    GlobalDataBase::UpdateData("inflowKnudsenNumber", &inflowKnudsenNumber, PHDOUBLE, 1);

    int nZones = PHMPI::GetNumberofGlobalZones();

    RDouble timeStep = 1.0e10;
    for(int iZone=0; iZone< nZones;iZone++)
    {
        //int iZone = ZoneNo[nZone];

        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);
        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));

        int nTotalCell = grid->GetNTotalCell();
        RDouble * cellVolume = grid->GetCellVolume();

        RDouble maxCellLength = 0.0;
        RDouble minCellLength = 1.0e10;
        RDouble aveCellLength = 0.0;

        for(int iCell = 0; iCell < nTotalCell; iCell++)
        {
            RDouble cellLength = pow(cellVolume[iCell],1.0/3.0);
            if(maxCellLength < cellLength)
            {
                maxCellLength = cellLength; 
            }else if(minCellLength > cellLength)
            {
                minCellLength = cellLength; 
            }
            aveCellLength = aveCellLength + cellLength;
        }
        aveCellLength = aveCellLength/nTotalCell;
        RDouble maxVelocity = inflowVelocity+1.0*sqrt(2.0*gasConstant*inflowTemperature);
        timeStep = min(timeStep, static_cast <RDouble> (minCellLength/3.0/maxVelocity));
    }

    //RDouble ftmp_min;
    //MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Allreduce(&ftmp, &ftmp_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD);

    timeStep = 0.95*min(timeStep,inflowMeanCollisionTime);

    oss << "Time step:"<<timeStep<<"\n"<<"inflowMeanCollisionTime:"<<inflowMeanCollisionTime<<endl;
    WriteLogFile(oss.str());

    this->timeStep = timeStep;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 函数头部说明：                                                             
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:   例如：ReadSpecies()                                             
// AUTHOR(S):  Li Jin                                                         
// PURPOSE:    读入DSMC计算所需组分参数信息                                    
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DSMC :: ReadSpecies()
{
    string filename = GlobalDataBase::GetStrParaFromDB("SpeciesData");
    fstream file;
    file.open(filename.c_str(),ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(filename);
    }

    nSpecies = GlobalDataBase::GetIntParaFromDB("nSpecies");

    species = new Species[nSpecies];

    string line, word;
    string seperator  = " =\t\r\n#$,;";

    do
    {
        ReadNewLine(file, line);
        if((line[0]) != '#') break;
    }   while(true); 

    for(int n = 0; n < nSpecies; n++)
    {
        line = FindNextWord(line,word,seperator);
        species[n].id = word;

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].percent,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].molwt,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].mass,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< int >(species[n].rotdof,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].rotrel,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< int >(species[n].vibdof,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].vibrel,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].vibtemp,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].specwt,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].charge,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].diam,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].omega,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].tref,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].alpha,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].rotc1,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].rotc2,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].vibc1,word,std::dec);

        line = FindNextWord(line,word,seperator);
        from_string< RDouble >(species[n].vibc2,word,std::dec);

        ReadNewLine(file, line);
    }

    file.close();

}
void DSMC::InitData()
{
    RDouble vstream[3];
    RDouble xcell[4],ycell[4],zcell[4];

    RDouble inflowVelocity = GlobalDataBase::GetDoubleParaFromDB("inflowVelocity");
    RDouble attackd    = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    RDouble attack   = attackd   * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    vstream[0] = inflowVelocity*cos(sideslip)*cos(attack);
    vstream[1] = inflowVelocity*cos(sideslip)*sin(attack);
    vstream[2] = inflowVelocity*sin(sideslip);


    int ispecies,jspecies;
    RDouble CCG1,CCG2;
    RDouble inflowTemperature = GlobalDataBase::GetDoubleParaFromDB("inflowTemperature");
    int nZones = PHMPI::GetNumberofGlobalZones();

    for(ispecies = 0; ispecies < nSpecies; ++ ispecies)
    {
        for(jspecies = 0; ispecies < nSpecies; ++ ispecies)
        {
            for(int izone = 0; izone < nZones; ++izone)
            {
                int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(izone);
                UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));
                int nTotalCell = grid->GetNTotalCell();  

                RDouble ** CCG = collision[izone]->GetCCG();

                for(int icell = 0; icell < nTotalCell; ++ icell)
                {
                    CCG1 = SPM[0][0][0]*300.*sqrt(inflowTemperature/300.);
                    CCG2 = myrand();
                    CCG[0][icell] = CCG1;
                    CCG[1][icell] = CCG2;
                }
            }
        }
    }


    int initialtype   = GlobalDataBase::GetIntParaFromDB("initialtype");
    nweight = GlobalDataBase::GetIntParaFromDB("nweight");
    RDouble inflowNumberDensity = GlobalDataBase::GetDoubleParaFromDB("inflowNumberDensity"); 

    Species * species = this->species;

    int nMolecule;
    if(initialtype > 0) 
    {
        for(int isp = 0; isp < nSpecies; ++isp)
        {
            if(species[isp].percent <= 1.e-20) continue; 
            
            RDouble remainder = 0.0;
            RDouble mostProbableMolecularThermalSpeed = sqrt(2.0*BOLTZ*inflowTemperature/species[isp].mass);

            for(int izone = 0; izone < nZones; ++izone)
            {
                int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(izone);
                UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));

                int nTotalCell = grid->GetNTotalCell();

                RDouble * cellVolume = grid->GetCellVolume();
                RDouble * cellWeight = dsmcgrid[izone]->GetCellWeight();
                int * cell2node = grid->GetCell2Node(); 

                RDouble * x = grid->GetX();
                RDouble * y = grid->GetY();
                RDouble * z = grid->GetZ();

                int * moleculeSpecies = molecule[izone]->GetMoleculeSpecies();
                int * moleculeLocation = molecule[izone]->GetMoleculeLocation();
                RDouble ** moleculePosition = molecule[izone]->GetMoleculePosition();
                RDouble ** moleculeVelocity = molecule[izone]->GetMoleculeVelocity();
                RDouble * moleculeRotationalEnergy = molecule[izone]->GetMoleculeRotationalEnergy();
                RDouble * moleculeVibrationalEnergy = molecule[izone]->GetMoleculeVibrationalEnergy();
                int * moleculeVibrationLevel = molecule[izone]->GetMoleculeVibrationLevel();

                RDouble Fnum = FNUMI[izone];

                if(isp == 0)
                {
                    nMolecule = 0;
                } 
                else
                {
                    nMolecule = dsmcgrid[izone]->GetnMolecule();
                }

                for(int icell = 0; icell < nTotalCell; ++icell)
                {
                    RDouble weight;
                    if(nweight == 0)
                    {
                        weight = 1.0;
                    } 
                    else
                    {
                        weight = cellWeight[icell];
                    }

                    RDouble numberofMolecule = inflowNumberDensity*cellVolume[icell]*species[isp].percent/(Fnum*weight) + remainder;

                    int integerNumberofMolecule;
                    if(izone == nZones-1 && icell == nTotalCell-1)
                    {
                        integerNumberofMolecule = int(numberofMolecule);
                    } 
                    else
                    {
                        integerNumberofMolecule = int(numberofMolecule);
                        remainder = numberofMolecule-integerNumberofMolecule;
                    }

                    if(integerNumberofMolecule > 0)
                    {
                        for(int iMolecule = 0; iMolecule < integerNumberofMolecule; ++ iMolecule)
                        {
                            if(nMolecule < MaxMoleculeNumber)
                            {
                                /*nMolecule += 1;*/

                                moleculeSpecies[nMolecule] = isp;
                                moleculeLocation[nMolecule] = icell;

                                for(int inode = 0; inode < 4; ++inode)
                                {
                                    xcell[inode] = x[cell2node[icell * 4 + inode]];
                                    ycell[inode] = y[cell2node[icell * 4 + inode]];
                                    zcell[inode] = z[cell2node[icell * 4 + inode]];
                                }

                                RDouble ss1 = myrand();
                                RDouble ss2 = sqrt(myrand());
                                RDouble ss3 = pow(myrand(),1.0/3.0);

                                RDouble xd = xcell[0] + ss1*(xcell[1] - xcell[0]);
                                RDouble yd = ycell[0] + ss1*(ycell[1] - ycell[0]);
                                RDouble zd = zcell[0] + ss1*(zcell[1] - zcell[0]);

                                RDouble xe = xcell[2] + ss2*(xd - xcell[2]);
                                RDouble ye = ycell[2] + ss2*(yd - ycell[2]);
                                RDouble ze = zcell[2] + ss2*(zd - zcell[2]);

                                moleculePosition[nMolecule][0] = xcell[3] +  ss3*(xe - xcell[3]);
                                moleculePosition[nMolecule][1] = ycell[3] +  ss3*(ye - ycell[3]);
                                moleculePosition[nMolecule][2] = zcell[3] +  ss3*(ze - zcell[3]);

                                RDouble vn = mostProbableMolecularThermalSpeed * sqrt(-log(myrand()));
                                RDouble vr = mostProbableMolecularThermalSpeed * sqrt(-log(myrand()));

                                RDouble theta1 = 2*PI*myrand();
                                RDouble theta2 = 2*PI*myrand();

                                moleculeVelocity[nMolecule][0] = vstream[0] + vn*cos(theta1);
                                moleculeVelocity[nMolecule][1] = vstream[1] + vr*cos(theta2);
                                moleculeVelocity[nMolecule][2] = vstream[2] + vr*sin(theta2);

                                moleculeRotationalEnergy[nMolecule] = molecule[izone]->CalMoleculeRotationalEnergy(species[isp].rotdof,inflowTemperature);

                                moleculeVibrationLevel[nMolecule] = static_cast <int> (-log(myrand())*inflowTemperature/species[isp].vibtemp);

                                moleculeVibrationalEnergy[nMolecule] = moleculeVibrationLevel[nMolecule]*BOLTZ*species[isp].vibtemp;
                                //molecule[izone]->CalMoleculeVibrationalEnergy(species[isp].vibdof,inflowTemperature,species[isp].vibtemp);
                                
                                nMolecule += 1;
                            }


                        
                        }
                    }
                }
                
                dsmcgrid[izone]->SetnMolecule(nMolecule);
            }
        }
    
    }

    //for(int izone = 0; izone < nZones; ++izone)
    //{
    //    cout<<"izone="<<izone<<", nMolecule="<<dsmcgrid[izone]->GetnMolecule()<<endl;
    //}



}

void DSMC::MoleculeIndex()
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    for(int izone = 0; izone < nZones; ++izone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(izone);
        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));
        int nTotalCell = grid->GetNTotalCell();
        
        int * numberOfMoleculeInCells = index[izone]->GetNumberOfMoleculeInCells();

        for(int icell = 0; icell < nTotalCell; ++ icell)
        {
            numberOfMoleculeInCells[icell] = 0;
        }
    }

    for(int izone = 0; izone < nZones; ++izone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(izone);
        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));
        int nTotalCell = grid->GetNTotalCell();
        
        int nMolecule = dsmcgrid[izone]->GetnMolecule();
        int * moleculeLocation = molecule[izone]->GetMoleculeLocation();

        int * startAddressOfMoleculeInCells = index[izone]->GetStartAddressOfMoleculeInCells();
        int * numberOfMoleculeInCells = index[izone]->GetNumberOfMoleculeInCells();
        int * moleculeNoInOrderOfCells = index[izone]->GetMoleculeNoInOrderOfCells();

        for(int imolecule = 0; imolecule < nMolecule; ++ imolecule)
        {
            int icell = moleculeLocation[imolecule];
            numberOfMoleculeInCells[icell] += 1;
        }

        int M = 0;
        for(int icell = 0; icell < nTotalCell; ++ icell)
        {
            startAddressOfMoleculeInCells[icell] = M;
            M = M + numberOfMoleculeInCells[icell];
        }

        for(int icell = 0; icell < nTotalCell; ++ icell)
        {
            numberOfMoleculeInCells[icell] = 0;
        }

        for(int imolecule = 0; imolecule < nMolecule; ++ imolecule)
        {
            int icell = moleculeLocation[imolecule];
            numberOfMoleculeInCells[icell] += 1;
            moleculeNoInOrderOfCells[imolecule] = startAddressOfMoleculeInCells[icell] + numberOfMoleculeInCells[icell];
        }
    }





}

void DSMC::InitSample()
{
    NSMP = 0;
    TIMI = TIME;

    int nZones = PHMPI::GetNumberofGlobalZones();
    nSpecies = GlobalDataBase::GetIntParaFromDB("nSpecies");
    MMVM = GlobalDataBase::GetIntParaFromDB("MMVM");

    for(int izone = 0; izone < nZones; ++izone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(izone);
        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));
        int nTotalCell = grid->GetNTotalCell();

        RDouble ** CSR = sample[izone]->GetCSR();
        RDouble *** CSV = sample[izone]->GetCSV();
        RDouble *** CS = sample[izone]->GetCS();

        for(int ispecies = 0; ispecies < nSpecies; ++ ispecies)
        {
            for(int icell = 0; icell < nTotalCell; ++ icell)
            {
                CSR[icell][ispecies] = 0.0;

                if(MMVM > 0) 
                {
                    //for(int K = 0; K < ISPV[ispecies]; ++ K)
                    for(int K = 0; K < 1; ++ K)
                    {
                        CSV[K][icell][ispecies] = 0.0;
                    }
                }

                for(int K = 0; K <=7; ++ K)
                {
                    CS[K][icell][ispecies] = 0.0;
                }

            }
        }
    }

    for(int izone = 0; izone < nZones; ++izone)
    {
        RDouble *** CSS1 = sample[izone]->GetCSS1();
        RDouble *** CSS2 = sample[izone]->GetCSS2();

        int nWallFace = dsmcgrid[izone]->GetnWallFace();

        if(nWallFace > 0) 
        {
            for(int ispecies = 0; ispecies < nSpecies; ++ ispecies)
            {
                for(int iface = 0; iface < nWallFace; ++ iface)
                {
                    for(int K = 0; K< 8; ++ K)
                    {
                        CSS1[K][iface][ispecies] = 0.0;
                        CSS2[K][iface][ispecies] = 0.0;
                    }
                }
            }
        }
    }   

}



void DSMC::MoleculeCollision()
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    int nweight = Getnweight();
    RDouble weight,ASEL;
    int NSEL,collision_pair_i,collision_pair_j;
    RDouble RML,RMM,m1plusm2;

    int MMVM = GetMMVM();

    for(int izone = 0; izone < nZones; ++izone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(izone);
        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));
        int nTotalCell = grid->GetNTotalCell();

        RDouble Fnum = FNUMI[izone];
        RDouble dtm = DTMI[izone];

        int * numberOfMoleculeInCells = index[izone]->GetNumberOfMoleculeInCells();

        //int * face_number_of_each_cell = grid->GetFaceNumberOfEachCell();
        RDouble * cellWeight = dsmcgrid[izone]->GetCellWeight();
        RDouble * cellVolume = grid->GetCellVolume();
        RDouble ** CCG = collision[izone]->GetCCG();

        for(int icell = 0; icell < nTotalCell; ++ icell)
        {
            int IC2 = numberOfMoleculeInCells[icell];

            if(IC2 > 1) 
            {
                if(nweight == 0)
                {
                    weight = 1.0;
                } 
                else
                {
                    weight = cellWeight[icell];
                }
                
                ASEL = 0.5*IC2*(IC2-1)*Fnum*weight*CCG[0][icell]*dtm/cellVolume[icell] + CCG[1][icell];
                NSEL = static_cast<int>(ASEL);
                CCG[1][icell] = ASEL - NSEL;

                RDouble ** CCG = collision[izone]->GetCCG();
                RDouble CVR = collision[izone]->GetCVR();
                vector<RDouble> VCM = collision[izone]->GetVCM();   
                int * moleculeSpecies = molecule[izone]->GetMoleculeSpecies();
                RDouble ** moleculeVelocity = molecule[izone]->GetMoleculeVelocity();
                int * moleculeVibrationLevel = molecule[izone]->GetMoleculeVibrationLevel();

                collision_pair_i = collision_pair_j = -1;
                if(NSEL > 0)
                {
                    for(int iselect = 0; iselect < NSEL; ++ iselect)
                    {
                        SelectMolecule(izone,icell,collision_pair_i,collision_pair_j);

                        if(CVR > CCG[0][icell])
                            CCG[0][icell] = CVR;
                        
                        if(CVR/CCG[0][icell] > myrand())
                        {
                            int isp = moleculeSpecies[collision_pair_i];
                            int jsp = moleculeSpecies[collision_pair_j];
                            if(isp<0 || jsp<0) 
                            {
                                cout<<isp<<" "<<jsp<<endl;
                                cin.get();
                            }
                            RML = SPM[4][isp][jsp]/species[jsp].mass;
                            RMM = SPM[4][isp][jsp]/species[isp].mass;
                            m1plusm2 = species[isp].mass + species[jsp].mass;

                            VCM[0] = RML*moleculeVelocity[0][collision_pair_i] + RMM*moleculeVelocity[0][collision_pair_j];
                            VCM[1] = RML*moleculeVelocity[1][collision_pair_i] + RMM*moleculeVelocity[1][collision_pair_j];
                            VCM[2] = RML*moleculeVelocity[2][collision_pair_i] + RMM*moleculeVelocity[2][collision_pair_j];


                            preCollisionEtran = 0.5*SPM[4][isp][jsp]*collision[izone]->GetVRR();
                            preCollisionErot = 0.0;
                            preCollisionEvib = 0.0;

                            if(species[isp].rotdof > 0) 
                                preCollisionErot += molecule[izone]->GetMoleculeRotationalEnergy()[collision_pair_i];
                            if(species[jsp].rotdof > 0) 
                                preCollisionErot += molecule[izone]->GetMoleculeRotationalEnergy()[collision_pair_j];

                            if(MMVM > 0) 
                            {
                                if(species[isp].vibdof > 0) 
                                    preCollisionEvib += moleculeVibrationLevel[collision_pair_i]*BOLTZ*species[isp].vibtemp;
                                if(species[jsp].vibdof > 0) 
                                    preCollisionEvib += moleculeVibrationLevel[collision_pair_j]*BOLTZ*species[jsp].vibtemp;
                            }
                            preCollisionEtotal = preCollisionEtran + preCollisionErot + preCollisionEvib;
                            preCollisionEint =  preCollisionErot + preCollisionEvib;

                            postCollisionEtotal = preCollisionEtotal;

                            int reactFlag = 0;
                            if(reactFlag == 1)
                            {
                            
                            }
                            else if(reactFlag == 0)
                            {
                                NonReactingEDisposal(izone,collision_pair_i,collision_pair_j);

                                TwobodyScattering(izone,collision_pair_i,collision_pair_j,m1plusm2,reactFlag);
                            
                            
                            }

                        }
                    }


                
                }
            }
        
        }
    }



}
void DSMC::TwobodyScattering(int izone, int collision_pair_i, int collision_pair_j, RDouble m1plusm2, int & reactFlag)
{
    RDouble VCM[3],VRC[3],VRCP[3];

    int isp = molecule[izone]->GetMoleculeSpecies()[collision_pair_i];
    int jsp = molecule[izone]->GetMoleculeSpecies()[collision_pair_j];

    RDouble RML = SPM[5][isp][jsp]/species[jsp].mass;
    RDouble RMM = SPM[5][isp][jsp]/species[isp].mass;

    RDouble m3plusm4 = species[isp].mass + species[jsp].mass;

    RDouble wRatio = sqrt(m1plusm2/m3plusm4);

    for(int K = 0; K < 3; ++K)
    {
        VCM[K] = collision[izone]->GetVCM()[K];
        VRC[K] = collision[izone]->GetVRC()[K];
    }

    if(postCollisionEtran < 0)
        postCollisionEtran = 1.0e-26;

    RDouble VR = sqrt(2.0*postCollisionEtran/SPM[5][isp][jsp]);

    RDouble A,B,C,D;
    if(abs(SPM[4][isp][jsp]-1.0) < 1.0e-3)
    {
        B = 2*myrand()-1.0;
        A = sqrt(1-B*B);
        VRCP[0] = B*VR;
        C = 2.0*PI*myrand();
        VRCP[1] = A*cos(C)*VR;
        VRCP[2] = A*sin(C)*VR;
    }
    else
    {
        B = 2.0*(pow(myrand(),SPM[4][isp][jsp])) - 1.0;
        A = sqrt(1-B*B);
        C = 2.0*PI*myrand();
        D = sqrt(VRC[1]*VRC[1] + VRC[2]*VRC[2]);
        if(D > 1.0e-6)
        {
            VRCP[0] = B*VRC[0] + A*sin(C)*D;
            VRCP[1] = B*VRC[1] + A*(VR*VRC[2]*cos(C) - VRC[0]*VRC[1]*sin(C))/D;
            VRCP[2] = B*VRC[2] + A*(VR*VRC[1]*cos(C) + VRC[0]*VRC[2]*sin(C))/D;
        }
        else
        {
            VRCP[0] = B*VRC[0];
            VRCP[1] = A*cos(C)*VRC[0];
            VRCP[2] = A*sin(C)*VRC[0];
        }
    
    }
    for(int K = 0; K < 3; ++ K)
    {
        molecule[izone]->GetMoleculeVelocity()[K][collision_pair_i] = wRatio*VCM[K] + VRCP[K]*RMM;
        molecule[izone]->GetMoleculeVelocity()[K][collision_pair_j] = wRatio*VCM[K] - VRCP[K]*RML;
    
    }
}
void DSMC::NonReactingEDisposal(int izone, int collision_pair_i, int collision_pair_j)
{
    RDouble E_Dispose = preCollisionEtran;

    postCollisionErot = 0.0;
    postCollisionEvib = 0.0;
    postCollisionEint = 0.0;

    Species * species = this->species;

    int collision_pair_k;

    for(int ksp = 1; ksp < 3; ++ ksp)
    {
        if(ksp == 1)
            collision_pair_k = collision_pair_i;
        else if(ksp == 2)
            collision_pair_k = collision_pair_j;

        int KS = molecule[izone]->GetMoleculeSpecies()[collision_pair_i];

        RDouble rotrel = species[KS].rotrel;
        RDouble rotdof = species[KS].rotdof;
        RDouble fraction ;
        if(rotdof > 0)
        {
            if(rotrel > myrand()) 
            {
                E_Dispose += molecule[izone]->GetMoleculeRotationalEnergy()[collision_pair_k];
                if(rotdof == 2)
                    fraction = 1.0-pow(myrand(),1.0/(2.5-species[KS].omega));
                else 
                    Larsen_Borgenakke(0.5*rotdof-1.0,1.5-species[KS].omega,fraction);
            
            molecule[izone]->GetMoleculeRotationalEnergy()[collision_pair_k] = fraction*E_Dispose;
            E_Dispose -= molecule[izone]->GetMoleculeRotationalEnergy()[collision_pair_k];
            }
        }
        postCollisionErot += molecule[izone]->GetMoleculeRotationalEnergy()[collision_pair_k];

        int MMVM = GetMMVM();
        RDouble vibref = species[KS].vibrel;
        RDouble vibdof = species[KS].vibdof;
        if(MMVM > 0) 
        {
            if(vibref > myrand())
            {
                if(vibdof == 2) 
                {
                    RDouble Evib = molecule[izone]->GetMoleculeVibrationLevel()[collision_pair_k]*BOLTZ*species[KS].vibtemp;
                    E_Dispose += Evib;
                    int MaxLev = static_cast<int> (E_Dispose/(BOLTZ*species[KS].vibtemp));
                    RDouble prob = 0.0;
                    while(prob < myrand())
                    {
                        int IV = static_cast<int> (myrand()*(MaxLev + 0.99999));
                        molecule[izone]->GetMoleculeVibrationLevel()[collision_pair_k] = IV;
                        Evib = IV*BOLTZ*species[KS].vibtemp;
                        if(Evib < E_Dispose)
                            prob = pow(1-Evib/E_Dispose,1.5-species[KS].omega);

                    }
                    E_Dispose -= Evib;
                }
            }
            postCollisionEvib += molecule[izone]->GetMoleculeVibrationLevel()[collision_pair_k]*BOLTZ*species[KS].vibtemp;

        }
    }

    postCollisionEint = postCollisionErot + postCollisionEvib;
    postCollisionEtran = E_Dispose;
    postCollisionEtotal = postCollisionEtran + postCollisionEint;

}

void DSMC::Larsen_Borgenakke(RDouble xma,RDouble xmb,RDouble & fraction)
{
    RDouble P;
    while(1)
    {
        fraction = myrand();
        if(xma < 1.0e-6 || xmb < 1.0e-6)
        {
            if(xma < 1.0e-6 && xmb < 1.0e-6) return;
            if(xma < 1.0e-6) P = pow((1.0-fraction),xmb);
            if(xmb < 1.0e-6) P = pow((1.0-fraction),xma);
        }
        else
            P = pow(((xma+xmb)*fraction/xma),xma)*pow((xma+xmb)*(1-fraction)/xmb,xmb);

        if(P > myrand()) 
            exit(0);
    }
}


void DSMC::SelectMolecule(int izone,int icell,int & collision_pair_i,int & collision_pair_j)
{
    int * startAddressOfMoleculeInCells = index[izone]->GetStartAddressOfMoleculeInCells();
    int * numberOfMoleculeInCells = index[izone]->GetNumberOfMoleculeInCells();
    int * moleculeNoInOrderOfCells = index[izone]->GetMoleculeNoInOrderOfCells();
    
    int k;
    if(numberOfMoleculeInCells[icell] == 2)
    {
        k = startAddressOfMoleculeInCells[icell];
        collision_pair_i = moleculeNoInOrderOfCells[k];
        k = 1 + startAddressOfMoleculeInCells[icell];
        collision_pair_j = moleculeNoInOrderOfCells[k];

    }
    else
    {
        k = static_cast<int>(myrand()*numberOfMoleculeInCells[icell] + startAddressOfMoleculeInCells[icell] + 0);
        collision_pair_i = moleculeNoInOrderOfCells[k];
        collision_pair_j = collision_pair_i;

        while(collision_pair_j == collision_pair_i)
        {
            k = static_cast<int>(myrand()*numberOfMoleculeInCells[icell] + startAddressOfMoleculeInCells[icell] + 0);
            collision_pair_j = moleculeNoInOrderOfCells[k];
        }
    }

    RDouble ** moleculeVelocity = molecule[izone]->GetMoleculeVelocity();

    vector<RDouble> VRC;
    VRC.resize(3);

    for(k = 0; k < 3; ++ k)
    {
        VRC[k] = moleculeVelocity[k][collision_pair_i] - moleculeVelocity[k][collision_pair_j];
    }
    collision[izone]->SetVRC(VRC);

    RDouble VRR = VRC[0]*VRC[0] + VRC[1]*VRC[1] + VRC[2]*VRC[2];
    RDouble VR = sqrt(VRR);

    int * moleculeSpecies = molecule[izone]->GetMoleculeSpecies();

    int isp = moleculeSpecies[collision_pair_i];
    int jsp = moleculeSpecies[collision_pair_j];

    if(isp<0 || jsp<0)
    {
        cout<<"lla--"<<isp<<" "<<jsp<<endl;
        cout<<collision_pair_i<<" "<<collision_pair_i<<endl;
        cin.get();
    }

    RDouble CVR = VR * SPM[0][isp][jsp]*pow(2.0*BOLTZ*SPM[1][isp][jsp]/(SPM[4][isp][jsp]*VRR), SPM[2][isp][jsp]-0.5)/ SPM[5][isp][jsp];

    collision[izone]->SetVRR(VRR);
    collision[izone]->SetVR(VR);
    collision[izone]->SetCVR(CVR);
    

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 函数头部说明：                                                              +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:   例如：RDouble gamma(RDouble X)                                    +
// AUTHOR(S):  Li Jin                                                          +
// PURPOSE:    计算伽马函数                                                    +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RDouble gamma(RDouble X)
{
    RDouble A,Y,GAM;
    A = 1.0;
    Y = X;
    if(Y < 1.0)
    {
        A = A/Y;
    }
    else
    {
        Y = Y-1;
        do 
        {
            A = A*Y;
            Y = Y-1;
        }while(Y>=1.0);
    }
    GAM = A*(1.-0.5748646*Y+0.9512363*pow(Y,2)-0.6998588*pow(Y,3)+ 0.4245549*pow(Y,4)-0.1010678*pow(Y,5));
    return GAM;
}

RDouble gamma(int X)
{
    RDouble A,Y,GAM;
    A = 1.0;
    Y = X;
    if(Y < 1.0)
    {
        A = A/Y;
    }
    else
    {
        Y = Y-1;
        do 
        {
            A = A*Y;
            Y = Y-1;
        }while(Y>=1.0);
    }
    GAM = A*(1.-0.5748646*Y+0.9512363*pow(Y,2)-0.6998588*pow(Y,3)+ 0.4245549*pow(Y,4)-0.1010678*pow(Y,5));
    return GAM;
}





void DSMC::ZoneToCPU()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();
    int nCPU = GetNumberOfProcessor();

    ZonePerCPU = nZones/nCPU;
    ZoneLeft = nZones%nCPU;

    if(ZoneLeft==0)
    {
        ZonePerCPU = ZonePerCPU - 1;
        ZoneLeft  = nCPU;
    }
    int ZonePerCPU1,ZoneStart;
    if(myid < ZoneLeft)
    {
        ZonePerCPU1 = ZonePerCPU+1;
        ZoneStart = myid*ZonePerCPU1;
    }
    else
    {
        ZonePerCPU1 = ZonePerCPU;
        ZoneStart = ZoneLeft*(ZonePerCPU+1)+(myid-ZoneLeft)*ZonePerCPU1;
    }

    int ZoneEnd = ZoneStart+ZonePerCPU1;

    SetZoneStart(ZoneStart);
    SetZoneEnd(ZoneEnd);

    ZoneNo.resize (nZones);

    for(int i=0;i<nZones;i++)
    {
        ZoneNo[i]=i;
    }

}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 函数头部说明：                                                              +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:   例如：RDouble myrand()                                           +
// AUTHOR(S):  Li Jin                                                          +
// PURPOSE:    计算随机数                                                      +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
RDouble myrand()
{ 
  return rand()/RDouble(RAND_MAX);
}
}

