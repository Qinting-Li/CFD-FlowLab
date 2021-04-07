#include "SpecDiffHybSolver.h"
#include "SpecDiffHybGrid.h"
#include "ExplicitDifferenceBoundary.h"
#include "CompactDifferenceFirstDerivative.h"
#include "CompactDifferenceSecondDerivative.h"
#include "PoissonSolver.h"
#include "CorrectWallDivergence.h"
#include "RescaleField.h"
#include "Statistics.h"
#include "LES.h"
#include "p3dfft.h"
#include "TK_Time.h"
#include "TK_Exit.h"
#include "TK_Warning.h"
#include "PHMpi.h"
#include "PHHeader.h"
#include "GlobalDataBase.h"
#include <iostream>
#include <sstream>
#include <cstdlib>

namespace PHSPACE
{

    SpecDiffHybSolver::SpecDiffHybSolver()
    {
    }

    SpecDiffHybSolver::~SpecDiffHybSolver()
    {
        delete nonlinearTermCoef; nonlinearTermCoef = NULL;
        delete pressureGradientTermCoef; pressureGradientTermCoef = NULL;
        delete convectiveTerm; convectiveTerm = NULL;
        delete complexLaplacianVelocityTerm; complexLaplacianVelocityTerm = NULL;
        delete complexPressureGradientTerm; complexPressureGradientTerm = NULL;

        delete complexViscousTerm; complexViscousTerm = NULL;
        delete complexNonlinearTerm; complexNonlinearTerm = NULL;
        delete complexVelocityGradient; complexVelocityGradient = NULL;
        delete realVelocityGradient; realVelocityGradient = NULL;
        delete complexVelocitySquare; complexVelocitySquare = NULL;
        delete realVelocitySquare; realVelocitySquare = NULL;

        delete realnut; realnut = NULL;
        delete realnut_00; realnut_00 = NULL;
        delete realnut3d; realnut3d = NULL;
        delete complexU; complexU = NULL;
        delete complexV; complexV = NULL;
        delete complexW; complexW = NULL;
        delete complexP; complexP = NULL;
        delete complexDivergenceU; complexDivergenceU = NULL;

        delete GridData; GridData = NULL;
        delete DiffBoundaryData; DiffBoundaryData = NULL;
        delete CompactDifferenceFirstDerivativeData; CompactDifferenceFirstDerivativeData = NULL;
        delete CompactDifferenceSecondDerivativeData; CompactDifferenceSecondDerivativeData = NULL;
        delete RescaleFieldData; RescaleFieldData = NULL;
        delete PoissonSolverData; PoissonSolverData = NULL;
        delete CorrectWallDivergenceData; CorrectWallDivergenceData = NULL;
        delete StatisticsData; StatisticsData = NULL; 
        delete LESData; LESData = NULL;
    }

    void SpecDiffHybSolver::Run()
    {
        double last_time[3];
        TIME_SPACE::StartTime( last_time );
        using namespace PHMPI;
        int myid = GetCurrentProcessorID();

        Init();

        int maxsimustep = GlobalDataBase::GetIntParaFromDB("maxsimustep");
        int nstepstat = GlobalDataBase::GetIntParaFromDB("nstepstat");

        int istep0 = istep;
        do {
            istep = istep + 1;

            if( (myid == 0) && (istep%100 == 0) )
            {
                cout << setw(20) << istep;
                cout << setiosflags(ios::right);
                cout << setprecision(10);
                cout << setiosflags(ios::scientific);
                cout << setiosflags(ios::showpoint);
                cout << setw(20) << real((*complexU)(2, 2, 2)) << setw(20) << imag((*complexU)(2, 2, 2)) << "\n" << endl;
            }
            GetRHS();

            int iVortProjection = GlobalDataBase::GetIntParaFromDB("iVortProjection");
            if( iVortProjection == 0 ) 
            {
                Solve();
            }
            else
            {
                TK_Exit::UnexpectedVarValue( "iVortProjection", iVortProjection );
            }

            int instart = GlobalDataBase::GetIntParaFromDB("instart");
            int timeOrderOfPressureTerm = GlobalDataBase::GetIntParaFromDB("timeOrderOfPressureTerm");
            if( instart == 0 && ((istep - istep0) == timeOrderOfPressureTerm) )
            {
                SetTimeDataCoef();
            }
            else if( (IsConvectiveTermFailed || IsPressureGradientFailed) && ((istep - istep0) == timeOrderOfPressureTerm) )
            {
                SetTimeDataCoef();
                IsConvectiveTermFailed = false;
                IsPressureGradientFailed = false;
            }

            RescaleVelocityField(istep);

            int nstepstat = GlobalDataBase::GetIntParaFromDB("nstepstat");
            if(istep >= nstepstat)
            {
                GetStatistics();
            }

            int nwrest = GlobalDataBase::GetIntParaFromDB("nwrest");
            if(istep % nwrest == 0)
            {
                OutputResidual(istep);

                OutputRestartData(istep);
            }

            IsInitialization = false;

            /* test SpecDiffHybSolver
            int istf = GridData -> iStartFourierIndex;
            int iedf = GridData -> iEndFourierIndex  ;
            int jstf = GridData -> jStartFourierIndex;
            int jedf = GridData -> jEndFourierIndex  ;
            int kstf = GridData -> kStartFourierIndex;
            int kedf = GridData -> kEndFourierIndex  ;
            stringstream strmMyid;
            strmMyid << myid;
            string test = "./test" +  strmMyid.str() + ".dat";;
            fstream filetest;
            filetest.open(test.c_str(), ios_base::out|ios_base::app);
            filetest << "complexUVWP = \n";
            for(int ix = kstf; ix <= kedf; ++ix)
            {
                for(int iy = jstf; iy <= jedf; ++iy)
                {
                    for(int iz = istf; iz <= iedf; ++iz)
                    {
                        filetest << setw(20) << ix
                                 << setw(20) << iy
                                 << setw(20) << iz
                                 << setw(20) << (*complexU)(iz, iy, ix)
                                 << setw(20) << (*complexV)(iz, iy, ix)
                                 << setw(20) << (*complexW)(iz, iy, ix)
                                 << setw(20) << (*complexP)(iz, iy, ix)
                                 << "\n";
                    }
                }
            }
            filetest.close();
            filetest.clear();
            test SpecDiffHybSolver*/
            
        } while(istep < maxsimustep);

    }

    void SpecDiffHybSolver::Init()
    {

        GridData = new SpecDiffHybGrid();
        GridData -> InitGridData();

        InitTimeDependentData();

        DiffBoundaryData = new ExplicitDifferenceBoundary();
        DiffBoundaryData -> GetSchemeCoefDiffBoundaryData();

        CompactDifferenceFirstDerivativeData = new CompactDifferenceFirstDerivative();
        CompactDifferenceFirstDerivativeData -> InitCompactDifferenceFirstDerivativeData(GridData, DiffBoundaryData);

        CompactDifferenceSecondDerivativeData = new CompactDifferenceSecondDerivative();
        CompactDifferenceSecondDerivativeData -> InitCompactDifferenceSecondDerivativeData(GridData);		
        
        InitGlobal();

        CompactDifferenceFirstDerivativeData -> GetDetaDz(GridData);

        RescaleFieldData = new RescaleField();
        RescaleFieldData -> InitRescaleData( GridData );		

        InitFieldData();        

        PoissonSolverData = new PoissonSolver();
        PoissonSolverData -> InitPoissonSolverData( GridData, DiffBoundaryData, CompactDifferenceSecondDerivativeData, realnut, realnut_00 );

        int instart = GlobalDataBase::GetIntParaFromDB("instart");
        if(instart != 0)
        {
            InputOldTimeData();
        }

        CorrectWallDivergenceData = new CorrectWallDivergence();
        CorrectWallDivergenceData -> InitCorrectWallDivergenceData( GridData, CompactDifferenceFirstDerivativeData, CompactDifferenceSecondDerivativeData, PoissonSolverData, realnut );

        if ( instart == 0 )
        {
            VanishingDIVproj();
        }

        StatisticsData = new Statistics();
        StatisticsData -> InitStatisticsData(GridData);

        int ivisc = GlobalDataBase::GetIntParaFromDB("ivisc");
        if(ivisc > 1)
        {
            LESData = new LES();
            LESData -> InitLESData(GridData);
        }

        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        int procIDWithKx0Ky0 = GridData -> procIDWithKx0Ky0;
        if(myid == procIDWithKx0Ky0)
        {
            OutputFieldName();
        }
    }

    void SpecDiffHybSolver::InitTimeDependentData()
    {
        InitTimeDependentFileName();

        AllocTimeDependentData();

        int instart = GlobalDataBase::GetIntParaFromDB("instart");
        if(instart == 0)
        {
            SetInitTimeDataCoef();
        }
        else if(instart == 1)
        {
            SetTimeDataCoef();
        }
        else
        {
            TK_Exit::UnexpectedVarValue( "instart", instart );
        }
    }

    void SpecDiffHybSolver::InitTimeDependentFileName()
    {
        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        stringstream strmMyid;
        strmMyid << myid;

        string outputdir = "./results";
        outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");

        nonlinearTermFile = outputdir + "/NonlinearTerm_" + strmMyid.str() + ".rsta"; 
        nonlinearTermBackupFile = outputdir + "/NonlinearTerm_" + strmMyid.str() + "_bk.rsta";

        pressureGradientTermFile = outputdir + "/PressureGradientTerm_" + strmMyid.str() + ".rsta";
        pressureGradientTermBackupFile = outputdir + "/PressureGradientTerm_" + strmMyid.str() + "_bk.rsta";

        laplacianVelocityTermFile = outputdir + "/LaplacianVelocity_" + strmMyid.str() + ".rsta";
        laplacianVelocityTermBackupFile = outputdir + "/LaplacianVelocity_" + strmMyid.str() + "_bk.rsta";

        timeInformationRestartFile = outputdir + "/Timeinfo.rsta";
    }

    void SpecDiffHybSolver::AllocTimeDependentData()
    {
        int ist = GridData -> iStartFourierIndex;
        int ied = GridData -> iEndFourierIndex  ;
        int jst = GridData -> jStartFourierIndex;
        int jed = GridData -> jEndFourierIndex  ;
        int kst = GridData -> kStartFourierIndex;
        int ked = GridData -> kEndFourierIndex;

        int timeOrderOfConvectiveTerm = GlobalDataBase::GetIntParaFromDB("timeOrderOfConvectiveTerm"); 
        Range I( -(timeOrderOfConvectiveTerm-1), 0 );
        nonlinearTermCoef = new RDouble1D( I, fortranArray );

        int timeOrderOfPressureTerm = GlobalDataBase::GetIntParaFromDB("timeOrderOfPressureTerm");
        Range J( -timeOrderOfPressureTerm, 0 );
        Range K( -timeOrderOfPressureTerm, -1);
        pressureGradientTermCoef = new RDouble1D( J, fortranArray );

        Range II( ist , ied );
        Range JJ( jst , jed );
        Range KK( kst , ked );
        Range NN( 1   , 3   );
        convectiveTerm = new Complex4D( II, JJ, KK, NN , fortranArray );
        complexLaplacianVelocityTerm = new Complex4D( II, JJ, KK, NN, fortranArray );
        complexPressureGradientTerm = new Complex5D( II, JJ, KK, NN, K, fortranArray );
    }

    void SpecDiffHybSolver::SetInitTimeDataCoef()
    {
        (*nonlinearTermCoef)( 0) = 1.0;
        (*nonlinearTermCoef)(-1) = 0.0;

        int timeOrderOfPressureTerm = GlobalDataBase::GetIntParaFromDB("timeOrderOfPressureTerm");
        for ( int i = -timeOrderOfPressureTerm ; i <= 0 ; ++i )
        {
            (*pressureGradientTermCoef)(i) = 0.0;
        }
    }

    void SpecDiffHybSolver::SetTimeDataCoef()
    {
        (*nonlinearTermCoef)( 0) = 1.5;
        (*nonlinearTermCoef)(-1) = -0.5;

        int timeOrderOfPressureTerm = GlobalDataBase::GetIntParaFromDB("timeOrderOfPressureTerm");
        if ( timeOrderOfPressureTerm == 1 )
        {
            (*pressureGradientTermCoef)( 0) = 0.0;
            (*pressureGradientTermCoef)(-1) = 1.0;
        }
        else if ( timeOrderOfPressureTerm == 2 )
        {
            (*pressureGradientTermCoef)( 0) = 0.0;
            (*pressureGradientTermCoef)(-1) = 2.0;
            (*pressureGradientTermCoef)(-2) = -1.0;
        }
        else
        {
            TK_Exit::UnexpectedVarValue( "timeOrderOfPressureTerm", timeOrderOfPressureTerm );
        }
    }

    void SpecDiffHybSolver::InitGlobal()
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex  ;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex  ;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex  ;
        Range IF(istf, iedf);
        Range JF(jstf, jedf);
        Range KF(kstf, kedf);

        int istp = GridData -> iStartPhysicalIndex;
        int iedp = GridData -> iEndPhysicalIndex  ;
        int jstp = GridData -> jStartPhysicalIndex;
        int jedp = GridData -> jEndPhysicalIndex  ;
        int kstp = GridData -> kStartPhysicalIndex;
        int kedp = GridData -> kEndPhysicalIndex  ;
        Range IP(istp, iedp);
        Range JP(jstp, jedp);
        Range KP(kstp, kedp);

        istep = 0;

        int instart = GlobalDataBase::GetIntParaFromDB("instart");
        if(instart == 0)
        {
            IsInitialization = true;
        }
        else
        {
            IsInitialization = false;
        }

        IsConvectiveTermFailed = false;
        IsPressureGradientFailed = false;

        screenOutFile = "screen.out";

        complexVelocityGradient = new Complex4D(IF, JF, KF, Range(1, 12), fortranArray);
        realVelocityGradient = new RDouble4D(IP, JP, KP, Range(1, 12), fortranArray);
        (*complexVelocityGradient)(IF, JF, KF, Range(1, 12)) = PHComplex(0.0, 0.0);
        (*realVelocityGradient)(IP, JP, KP, Range(1, 12)) = 0.0;

        complexVelocitySquare = new Complex4D(IF, JF, KF, Range(1, 12), fortranArray);
        realVelocitySquare = new RDouble4D(IP, JP, KP, Range(1, 12), fortranArray);
        (*complexVelocitySquare)(IF, JF, KF, Range(1, 12)) = PHComplex(0.0, 0.0);
        (*realVelocitySquare)(IP, JP, KP, Range(1, 12)) = 0.0; 

        complexViscousTerm = new Complex4D(IF, JF, KF, Range(1, 3), fortranArray);
        (*complexViscousTerm)(IF, JF, KF, Range(1, 3)) = PHComplex(0.0, 0.0);
        complexNonlinearTerm = new Complex4D(IF, JF, KF, Range(1, 3), fortranArray);
        (*complexNonlinearTerm)(IF, JF, KF, Range(1, 3)) = PHComplex(0.0, 0.0);
    }

    void SpecDiffHybSolver::InitFieldData()
    {
        InitFieldFileName();

        AllocFieldData();

        GetInitField();
    }

    void SpecDiffHybSolver::InitFieldFileName()
    {
        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        stringstream strmMyid;
        strmMyid << myid;

        string outputdir = "./results";
        outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");
;
        velocityRestartFile = outputdir + "/Velocity_" + strmMyid.str() + ".rsta";
        velocityBackupFile = outputdir + "/Velocity_" + strmMyid.str() + "_bk.rsta";
        pressureRestartFile = outputdir + "/Pressure_" + strmMyid.str() + ".rsta";
        pressureBackupFile = outputdir + "/Pressure_" + strmMyid.str() + "_bk.rsta";
        velocityNameFile = outputdir + "/Velocity.nam";
        fieldNameFile = outputdir + "/field.nam";

    }

    void SpecDiffHybSolver::AllocFieldData()
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex  ;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex  ;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex  ;

        int istp = GridData -> iStartPhysicalIndex;
        int jstp = GridData -> jStartPhysicalIndex;
        int kstp = GridData -> kStartPhysicalIndex;
        int iedp = GridData -> iEndPhysicalIndex;
        int jedp = GridData -> jEndPhysicalIndex;
        int kedp = GridData -> kEndPhysicalIndex;

        Range IF(istf, iedf);
        Range JF(jstf, jedf);
        Range KF(kstf, kedf);

        Range IP(istp, iedp);
        Range JP(jstp, jedp);
        Range KP(kstp, kedp);

        complexU           = new Complex3D(IF, JF, KF, fortranArray);
        complexV           = new Complex3D(IF, JF, KF, fortranArray);
        complexW           = new Complex3D(IF, JF, KF, fortranArray);
        complexP           = new Complex3D(IF, JF, KF, fortranArray);
        complexDivergenceU = new Complex3D(IF, JF, KF, fortranArray);
        (*complexU)(IF, JF, KF) = PHComplex(0.0, 0.0);
        (*complexV)(IF, JF, KF) = PHComplex(0.0, 0.0);
        (*complexW)(IF, JF, KF) = PHComplex(0.0, 0.0);
        (*complexP)(IF, JF, KF) = PHComplex(0.0, 0.0);
        (*complexDivergenceU)(IF, JF, KF) = PHComplex(0.0, 0.0);

        realnut = new RDouble1D(IF, fortranArray);
        realnut_00 = new RDouble1D(IF, fortranArray);
        realnut3d = new RDouble3D(IP, JP, KP, fortranArray);
        (*realnut)(IF) = 0.0;
        (*realnut_00)(IF) = 0.0;
        (*realnut3d)(IP, JP, KP) = 0.0;
    }

    void SpecDiffHybSolver::GetInitField()
    {
        using namespace PHMPI;
        int myid = GetCurrentProcessorID();

        if(myid == server)
        {
            WriteVelocityName();
        }

        int instart = GlobalDataBase::GetIntParaFromDB("instart");
        if( instart == 0 )
        {
            InitProfile();
        }
        else if( instart == 1 )
        {
            InitFromOldField();
        }
        else
        {
            TK_Exit::UnexpectedVarValue( "instart", instart );
        }

        OutputFieldRestartFile(0);
    }

    void SpecDiffHybSolver::WriteVelocityName()
    {
        fstream file;
        file.open(velocityNameFile.c_str(), ios_base::out);
        if ( !file )
        {
            TK_Exit::ExceptionExit("could not open velocityNameFile\n");
        }
        file << "u;velocity \n";
        file << "v \n";
        file << "w \n";
        file.close();
        file.clear();
    }

    void SpecDiffHybSolver::InputOldTimeData()
    {
        ReadConvectiveTerm();
        ReadTimeDependentInformation();
        if(IsConvectiveTermFailed || IsPressureGradientFailed)
        {
            SetInitTimeDataCoef();
        }
    }

    void SpecDiffHybSolver::ReadConvectiveTerm()
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex;
        Range I(istf, iedf);
        Range J(jstf, jedf);
        Range K(kstf, kedf);
        Range N(1, 3);

        fstream file;
        file.open(nonlinearTermFile.c_str(), ios_base::in|ios_base::binary);
        if(!file)
        {
            IsConvectiveTermFailed = true;
            (*convectiveTerm)(I, J, K, N) = PHComplex(0.0, 0.0);     
            TK_Warning::Warning("ReadConvectiveTerm: NonlinearTerm.rsta not exist. ConvectiveTerm will be set zero.");
        }
        else
        {
            int inttmp[5];
            int nin = 3 * (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);
            file.read(reinterpret_cast<char*>(inttmp), 5*sizeof(int));
            file.read(reinterpret_cast<char*>(&(*convectiveTerm)(istf, jstf, kstf, 1)), nin*sizeof(PHComplex));
        }
        file.close();
        file.clear();
    }

    void SpecDiffHybSolver::ReadTimeDependentInformation()
    {
        using namespace PHMPI;
        int myid = GetCurrentProcessorID();

        fstream file;        
        int procIDWithKx0Ky0 = GridData -> procIDWithKx0Ky0;
        if(myid == procIDWithKx0Ky0)
        {
            file.open(timeInformationRestartFile.c_str(), ios_base::in);
            if(IsInitialization || (!file))
            {
                istep = 0;
            }
            else
            {
                file >> istep;
            }
            file.close();
            file.clear();
        }
#ifdef PH_PARALLEL
        MPI_Bcast(&istep, 1, MPI_INTEGER, procIDWithKx0Ky0, MPI_COMM_WORLD); 
#endif
    }

    void SpecDiffHybSolver::InitProfile()
    {      
        int ist = GridData -> iStartFourierIndex;
        int ied = GridData -> iEndFourierIndex;
        int jst = GridData -> jStartFourierIndex;
        int jed = GridData -> jEndFourierIndex;
        int kst = GridData -> kStartFourierIndex;
        int ked = GridData -> kEndFourierIndex;
        Range I(ist, ied);
        Range J(jst, jed);
        Range K(kst, ked);
        RDouble1D *uAvg = new RDouble1D(I, fortranArray);
        Complex1D *cfu = new Complex1D(I, fortranArray);
        Complex1D *cfv = new Complex1D(I, fortranArray);
        Complex1D *cfw = new Complex1D(I, fortranArray);
        Complex1D *cfp = new Complex1D(I, fortranArray);

        const double pi = 2.0 * acos(0.0);
        const int nXZ = (jed - jst + 1) * (ked - kst + 1);

        (*complexP)(I, J, K) = PHComplex(0.0, 0.0);
        (*complexU)(I, J, K) = PHComplex(0.0, 0.0);
        (*complexV)(I, J, K) = PHComplex(0.0, 0.0);
        (*complexW)(I, J, K) = PHComplex(0.0, 0.0);

        RDouble1D *rZ = GridData -> realZ;
        for( int iz = ist; iz <= ied; ++iz )
        {
            (*uAvg)(iz) = 1.5 * (1.0 - (*rZ)(iz) * (*rZ)(iz));
        }
        double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");
        double uavgDivideByUtau = 2.960975610 + 2.439024390 * log(reynolds);
        (*uAvg)(I) = (*uAvg)(I) * uavgDivideByUtau;

        RDouble1D *realKX = GridData -> realKX;
        RDouble1D *realKY = GridData -> realKY;
        double kNyquistX = GridData -> realKNyquistX;
        double kNyquistY = GridData -> realKNyquistY;
        double rrrmax = -1.0e9;
        double ekx = 1.0;
        double eky = 1.0;
        double ek = 1.0;

        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        srand( (unsigned int)time(NULL) + (unsigned int)myid );
        for( int ix = kst; ix <= ked ;++ix )
        {
            for( int iy = jst; iy <= jed; ++iy )
            {
                for( int iz = ist; iz <= ied; ++iz)
                {
                    double kx = (*realKX)(ix);
                    double ky = (*realKY)(iy);
                    double k2  = kx * kx + ky * ky;

                    double rd = abs( (*rZ)(iz) - (*rZ)(ist) );
                    rd = min( rd, abs( (*rZ)(iz) - (*rZ)(ied) ) );
                    double rdplus = rd * reynolds;

                    // !ekx
                    double kt = abs(kx);
                    if( kt > TINY )
                    {
                        if( rdplus > 100.0)
                        {
                            GetAssumedSpectrum(kt, rd, ekx);
                        }
                        else if( rd > TINY )
                        {
                            GetAssumedSpectrumNearWall(kt, rd, ekx);
                        }
                        else
                        {
                            ekx = 1.0;
                        }
                    }
                    else
                    {
                        ekx = 1.0;
                    }

                    // !eky
                    kt = abs(ky);
                    if( kt > TINY )
                    {
                        if( rdplus > 100.0)
                        {
                            GetAssumedSpectrum(kt, rd, eky);
                        }
                        else if( rd > TINY )
                        {
                            GetAssumedSpectrumNearWall(kt, rd, eky);
                        }
                        else
                        {
                            eky = 1.0;
                        }
                    }
                    else
                    {
                        eky = 1.0;
                    }

                    ek = ekx * eky;
                    double rrr = ek;

                    if( rrr > rrrmax )
                    {
                        rrrmax = rrr;
                    }

                    if( (abs(kx) < kNyquistX) && (abs(ky) < kNyquistY) )
                    {
                        double theta = rand()/double(RAND_MAX);
                        theta = 2.0 * pi * theta;
                        double zz = rand()/double(RAND_MAX);
                        zz = 2.0 * zz - 1.0;
                        (*complexU)(iz, iy, ix) = rrr * PHComplex(cos(theta), sin(theta)) * zz;

                        theta = rand()/double(RAND_MAX);
                        theta = 2.0 * pi * theta;
                        zz = rand()/double(RAND_MAX);
                        zz = 2.0 * zz - 1.0;
                        (*complexV)(iz, iy, ix) = rrr * PHComplex(cos(theta), sin(theta)) * zz;

                        theta = rand()/double(RAND_MAX);
                        theta = 2.0 * pi * theta;
                        zz = rand()/double(RAND_MAX);
                        zz = 2.0 * zz - 1.0;
                        (*complexW)(iz, iy, ix) = rrr * PHComplex(cos(theta), sin(theta)) * zz;
                    }
                    else
                    {
                        continue;
                    }
                }// !for iz
            }// !for iy
        }// !for ix

        if( rrrmax > 1.0/nXZ )
        {
            (*complexU)(I, J, K) = (*complexU)(I, J, K) / rrrmax * 1.0 / static_cast<double>(nXZ);
            (*complexV)(I, J, K) = (*complexV)(I, J, K) / rrrmax * 1.0 / static_cast<double>(nXZ);
            (*complexW)(I, J, K) = (*complexW)(I, J, K) / rrrmax * 1.0 / static_cast<double>(nXZ);
        }
        
        // !filtering the fields:
        for(int ix = kst; ix <= ked; ++ix)
        {
            for(int iy = jst; iy <= jed; ++iy)
            {
                (*cfu)(I) = (*complexU)(I, iy, ix);
                (*cfv)(I) = (*complexV)(I, iy, ix);
                (*cfw)(I) = (*complexW)(I, iy, ix);
                (*cfp)(I) = (*complexP)(I, iy, ix);
                for(int iz = ist + 1; iz <= ied - 1; ++iz)
                {
                    (*complexU)(iz, iy, ix) = ( (*cfu)(iz - 1) + (*cfu)(iz + 1) ) / 6.0 + 2.0 * (*cfu)(iz) / 3.0;
                    (*complexV)(iz, iy, ix) = ( (*cfv)(iz - 1) + (*cfv)(iz + 1) ) / 6.0 + 2.0 * (*cfv)(iz) / 3.0;
                    (*complexW)(iz, iy, ix) = ( (*cfw)(iz - 1) + (*cfw)(iz + 1) ) / 6.0 + 2.0 * (*cfw)(iz) / 3.0;
                    (*complexP)(iz, iy, ix) = ( (*cfp)(iz - 1) + (*cfp)(iz + 1) ) / 6.0 + 2.0 * (*cfp)(iz) / 3.0;
                }
            }
        }
        
        //! clearing (0,0) mode:
        int x0Fall = GridData -> x0Fall;
        int y0Fall = GridData -> y0Fall;
        if( (kst == x0Fall) && (jst == y0Fall) )
        {
            (*complexU)(I, jst, kst) = PHComplex(0.0, 0.0);
            (*complexV)(I, jst, kst) = PHComplex(0.0, 0.0);
            (*complexW)(I, jst, kst) = PHComplex(0.0, 0.0);
            (*complexP)(I, jst, kst) = PHComplex(0.0, 0.0);
        }
        
        // !vanishing velocity on walls:
        (*complexU)(ist, J, K) = PHComplex(0.0, 0.0);
        (*complexV)(ist, J, K) = PHComplex(0.0, 0.0);
        (*complexW)(ist, J, K) = PHComplex(0.0, 0.0);

        (*complexU)(ied, J, K) = PHComplex(0.0, 0.0);
        (*complexV)(ied, J, K) = PHComplex(0.0, 0.0);
        (*complexW)(ied, J, K) = PHComplex(0.0, 0.0);

        // !cu(iz,jst,kst) = cmplx(uavg,0)
        if( (jst == y0Fall) && (kst == x0Fall) )
        {
            for( int iz = ist; iz <= ied; ++iz )
            {
                (*complexP)(iz, jst, kst) = PHComplex(0.0, 0.0);
                (*complexU)(iz, jst, kst) = PHComplex((*uAvg)(iz), 0.0);
            }
        }
        
        RescaleVelocityField(0);

        delete uAvg; uAvg = NULL;
        delete cfu; cfu = NULL;
        delete cfv; cfv = NULL;
        delete cfw; cfw = NULL;
        delete cfp; cfp = NULL;
    }

    void SpecDiffHybSolver::InitFromOldField()
    {
        ReadinVelocityRestartFile();

        int iVortProjection = GlobalDataBase::GetIntParaFromDB("iVortProjection");
        if( iVortProjection == 0 )
        {
            ReadinPressureRestartFile();
        }
        else
        {
            return;
        }
    }

    void SpecDiffHybSolver::OutputFieldRestartFile(int istep)
    {
        int nwrest = GlobalDataBase::GetIntParaFromDB("nwrest");
        int nTwrest = istep / nwrest;
        
        string velocityFileName = "Velocity.rsta";
        string pressureFileName = "Pressure.rsta";
        if( nTwrest % 2 == 0 )
        {
            velocityFileName = velocityBackupFile;
            pressureFileName = pressureBackupFile;
        }
        else
        {
            velocityFileName = velocityRestartFile;
            pressureFileName = pressureRestartFile;
        }

        WriteVelocityRestartFile(velocityFileName);
        int iVortProjection = GlobalDataBase::GetIntParaFromDB("iVortProjection");
        if(iVortProjection == 0)
        {
            WritePressureRestartFile(pressureFileName);
        }
    }

    void SpecDiffHybSolver::GetAssumedSpectrum( const double kt, const double rd, double &ek )
    {
      const double p0 = 2.0;
      const double cL = 6.78;
      const double beta = 5.2;
      const double ceta = 0.4;

      double kd = kt*rd;
      double kL = kt * 3.0;
      const double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");
      double keta = kt / (0.1 * reynolds);

      double ekInertia = pow(kd, -5.0/3.0);
      double ekL = pow(kL / sqrt(kL * kL + cL), 5.0/3.0 + p0);
      double ekEpsilon = exp( -beta * ( pow(pow(keta, 4.0) + pow(ceta, 4.0), 0.25) - ceta ) );

      ek = ekEpsilon * ekL * ekInertia;
    }

    void SpecDiffHybSolver::GetAssumedSpectrumNearWall( const double kt, const double rd, double &ek )
    {
        const double beta = 5.2;
        double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");
        double rdplus = rd * reynolds;
        double keta = kt / (0.1 * rdplus); 

        double ekEpsilon = exp( -beta * keta );
        ek = pow(0.5, ekEpsilon);
    }
    
    void SpecDiffHybSolver::ReadinVelocityRestartFile()
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex;

        fstream file;
        file.open(velocityRestartFile.c_str(),ios_base::in|ios_base::binary);
        if ( !file )
        {
            TK_Exit::FileOpenErrorExit(velocityRestartFile);
        }

        int inttmp[5];
        file.read(reinterpret_cast<char*>(inttmp), 5*sizeof(int));
        int nin = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);
        file.read(reinterpret_cast<char*>(&(*complexU)(istf, jstf, kstf)), nin*sizeof(PHComplex));
        file.read(reinterpret_cast<char*>(&(*complexV)(istf, jstf, kstf)), nin*sizeof(PHComplex));
        file.read(reinterpret_cast<char*>(&(*complexW)(istf, jstf, kstf)), nin*sizeof(PHComplex));
        file.close();
        file.clear();
    }

    void SpecDiffHybSolver::ReadinPressureRestartFile()
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex;

        fstream file;
        file.open(pressureRestartFile.c_str(),ios_base::in|ios_base::binary);
        if ( !file )
        {
            TK_Exit::FileOpenErrorExit(pressureRestartFile);
        }

        int inttmp[5];
        file.read(reinterpret_cast<char*>(inttmp), 5*sizeof(int));
        int nin = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);
        file.read(reinterpret_cast<char*>(&(*complexP)(istf, jstf, kstf)), nin*sizeof(PHComplex));
        file.close();
        file.clear();
    }

    void SpecDiffHybSolver::WriteVelocityRestartFile(string velocityFileName)
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex;

        fstream file;
        file.open(velocityFileName.c_str(), ios_base::out|ios_base::binary);
        if ( !file )
        {
            TK_Exit::FileOpenErrorExit(velocityFileName);
        }

        int inttmp = 1;
        file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
        inttmp = 3;
        Int1D *fourierIndexSize = GridData -> fourierIndexSize;
        file.write(reinterpret_cast<char*>(&(*fourierIndexSize)(1)), 3*sizeof(int));
        file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
        int nout = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);
        file.write(reinterpret_cast<char*>(&(*complexU)(istf, jstf, kstf)), nout*sizeof(PHComplex));
        file.write(reinterpret_cast<char*>(&(*complexV)(istf, jstf, kstf)), nout*sizeof(PHComplex));
        file.write(reinterpret_cast<char*>(&(*complexW)(istf, jstf, kstf)), nout*sizeof(PHComplex));
        file.close();
        file.clear();
    }

    void SpecDiffHybSolver::WritePressureRestartFile( string pressureFileName )
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex;

        fstream file;
        file.open(pressureFileName.c_str(), ios_base::out|ios_base::binary);
        if (!file)
        {
            TK_Exit::FileOpenErrorExit(pressureFileName);
        }

        int inttmp = 1;
        file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
        Int1D *fourierIndexSize = GridData -> fourierIndexSize;
        file.write(reinterpret_cast<char*>(&(*fourierIndexSize)(1)), 3*sizeof(int));
        file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
        int nout = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);
        file.write(reinterpret_cast<char*>(&(*complexP)(istf, jstf, kstf)), nout*sizeof(PHComplex));
        file.close();
        file.clear();
    }

    void SpecDiffHybSolver::GetRHS()
    {
        GetVelocityGradient();
        GetNonlinearTerm();

        int ivisc = GlobalDataBase::GetIntParaFromDB("ivisc");
        int iVortProjection = GlobalDataBase::GetIntParaFromDB("iVortProjection");
        if( (ivisc != 0) || (iVortProjection !=0 ) )
        {
            GetViscousTerm();
        }		
    }

    void SpecDiffHybSolver::Solve()
    {
        const PHComplex cZero = PHComplex(0.0, 0.0);
        const PHComplex cOne = PHComplex(1.0, 0.0);
        const PHComplex cOneI = PHComplex(0.0, 1.0);

        RDouble2D *coef_LHS = CompactDifferenceSecondDerivativeData -> lHSCoefMatrix;
        RDouble2D *coef_RHS = CompactDifferenceSecondDerivativeData -> rHSCoefMatrix;
        Int1D *lowerBound = CompactDifferenceSecondDerivativeData -> lowerBoundofRHSCoefMatrix;
        Int1D *upperBound = CompactDifferenceSecondDerivativeData -> upperBoundofRHSCoefMatrix;
        RDouble4D *coefVelocity = PoissonSolverData -> coefMatrixVelocityDirichlet;

        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex  ;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex  ;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex  ;
        int nSolve = iedf - istf;
        Range IF(istf, iedf);
        RDouble1D *realPHI = new RDouble1D(IF , fortranArray);
        double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");
        (*realPHI)(IF) = 2.0 * reynolds/( 1.0 + (*realnut)(IF) );

        Complex1D *complexULocal = new Complex1D(IF , fortranArray);
        Complex1D *complexVLocal = new Complex1D(IF , fortranArray);
        Complex1D *complexWLocal = new Complex1D(IF , fortranArray);
        Complex1D *complexPLocal = new Complex1D(IF , fortranArray);
        Complex1D *cR = new Complex1D(IF, fortranArray);
        RDouble2D *rC = new RDouble2D(IF, Range(-2, 2), fortranArray);
        RDouble1D *realKX = GridData -> realKX;
        RDouble1D *realKY = GridData -> realKY;
        Complex1D *complexKX = GridData -> complexKX;
        Complex1D *complexKY = GridData -> complexKY;
        int x0Fall = GridData -> x0Fall;
        int y0Fall = GridData -> y0Fall;
        double rdp_add = GlobalDataBase::GetDoubleParaFromDB("rdp_add");

        if( (kstf == x0Fall) && (jstf == y0Fall) )
        {
            int imode00fix = GlobalDataBase::GetIntParaFromDB("imode00fix");
            if( imode00fix > 0 )
            {
                (*realPHI)(IF) = 2.0 * reynolds/( 1.0 + (*realnut_00)(IF) );
            }

            int ix = kstf;
            int iy = jstf;
            double kxx = (*realKX)(ix);
            double kyy = (*realKY)(iy);
            double k2  = kxx * kxx + kyy * kyy;

            for( int iz = istf; iz <= iedf; ++iz )
            {
                (*complexULocal)(iz) = -(*complexNonlinearTerm)(iz, iy, ix, 1) - PHComplex(rdp_add, 0.0) + (*complexViscousTerm)(iz, iy, ix, 1);
            }
            (*complexP)(IF, iy, ix) = cZero;

            PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(istf+1, -2), &(*complexULocal)(istf), &(*lowerBound)(istf+1), &(*upperBound)(istf+1), cZero, cZero, &(*realPHI)(istf), 
                                                  &(*coefVelocity)(istf, -2, iy, ix), &(*cR)(istf) );

            (*complexULocal)(IF) = (*complexU)(IF, iy, ix) + (*complexULocal)(IF);

            (*complexU)(IF, iy ,ix) = (*complexULocal)(IF);
            (*complexV)(IF, iy ,ix) = cZero;
            (*complexW)(IF, iy ,ix) = cZero;
        }

        (*realPHI)(IF) = 2.0 * reynolds/( 1.0 + (*realnut)(IF) );

        RDouble1D *detaDz = GridData -> realDetaDz;

        int nBoundary_Neumann = DiffBoundaryData -> nBoundary_Neumann;

        RDouble4D *coefMatrixPressureNeumann = PoissonSolverData -> coefMatrixPressureNeumann;
        RDouble4D *boundaryMatrixPressureNeumann = PoissonSolverData -> boundaryMatrixPressureNeumann;

        Complex3D *complexDeltaU1 = CorrectWallDivergenceData -> complexDeltaU1;
        Complex3D *complexDeltaU2 = CorrectWallDivergenceData -> complexDeltaU2;
        Complex3D *complexDeltaV1 = CorrectWallDivergenceData -> complexDeltaV1;
        Complex3D *complexDeltaV2 = CorrectWallDivergenceData -> complexDeltaV2;
        Complex3D *complexDeltaW1 = CorrectWallDivergenceData -> complexDeltaW1;
        Complex3D *complexDeltaW2 = CorrectWallDivergenceData -> complexDeltaW2;
        Complex3D *complexDeltaP1 = CorrectWallDivergenceData -> complexDeltaP1;
        Complex3D *complexDeltaP2 = CorrectWallDivergenceData -> complexDeltaP2;
        Complex2D *complexDivDeltaU1LowerBound = CorrectWallDivergenceData -> complexDivDeltaU1LowerBound;
        Complex2D *complexDivDeltaU1UpperBound = CorrectWallDivergenceData -> complexDivDeltaU1UpperBound;
        Complex2D *complexDivDeltaU2LowerBound = CorrectWallDivergenceData -> complexDivDeltaU2LowerBound;
        Complex2D *complexDivDeltaU2UpperBound = CorrectWallDivergenceData -> complexDivDeltaU2UpperBound;

        Complex1D *dVardZ = new Complex1D(IF, fortranArray);
        Complex1D *omegaZ = new Complex1D(IF, fortranArray);
        Complex1D *deltaPressure = new Complex1D(IF, fortranArray);
        Complex1D *divergenceOfUStar = new Complex1D(IF, fortranArray);        
        for( int ix = kstf; ix <= kedf; ++ix )
        {
            for( int iy = jstf; iy <= jedf; ++iy )
            {
                double kxx = (*realKX)(ix);
                double kyy = (*realKY)(iy);
                double k2  = kxx * kxx + kyy * kyy;
                PHComplex ckx = (*complexKX)(ix);
                PHComplex cky = (*complexKY)(iy);

                if( k2 < TINY )
                {
                    continue;
                }

                // ! Nyquist wavenumbers:
                double realKNyquistY = GridData -> realKNyquistY;
                double realKNyquistX = GridData -> realKNyquistX;
                if( (abs(kyy) == abs(realKNyquistY)) || (abs(kxx) == abs(realKNyquistX)) )
                {
                    (*complexP)(IF, iy, ix) = cZero;
                    (*complexU)(IF, iy ,ix) = cZero;
                    (*complexV)(IF, iy, ix) = cZero;
                    (*complexW)(IF, iy, ix) = cZero;
                    continue;
                }

                (*complexPLocal)(IF) = (*complexP)(IF, iy ,ix);
                CompactDifferenceFirstDerivativeData -> GetDvarDz( istf, iedf, &(*detaDz)(istf), &(*complexPLocal)(istf), &(*dVardZ)(istf), &(*cR)(istf), &(*rC)(istf, -2) );

                for( int iz = istf; iz <= iedf; ++iz)
                {
                    (*complexULocal)(iz) = -(*complexNonlinearTerm)(iz, iy, ix, 1) - ckx * (*complexPLocal)(iz) + (*complexViscousTerm)(iz, iy ,ix, 1);
                    (*complexVLocal)(iz) = -(*complexNonlinearTerm)(iz, iy, ix, 2) - cky * (*complexPLocal)(iz) + (*complexViscousTerm)(iz, iy ,ix, 2);
                    (*complexWLocal)(iz) = -(*complexNonlinearTerm)(iz, iy, ix, 3) - (*dVardZ)(iz) + (*complexViscousTerm)(iz, iy ,ix, 3);
                }

                 PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(istf+1, -2), &(*complexULocal)(istf), &(*lowerBound)(istf+1), &(*upperBound)(istf+1), cZero, cZero, &(*realPHI)(istf), 
                                                       &(*coefVelocity)(istf, -2, iy, ix), &(*cR)(istf) );
                 PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(istf+1, -2), &(*complexVLocal)(istf), &(*lowerBound)(istf+1), &(*upperBound)(istf+1), cZero, cZero, &(*realPHI)(istf), 
                                                       &(*coefVelocity)(istf, -2, iy, ix), &(*cR)(istf) );
                 PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(istf+1, -2), &(*complexWLocal)(istf), &(*lowerBound)(istf+1), &(*upperBound)(istf+1), cZero, cZero, &(*realPHI)(istf), 
                                                       &(*coefVelocity)(istf, -2, iy, ix), &(*cR)(istf) );

                 (*complexULocal)(IF) = (*complexU)(IF, iy, ix) + (*complexULocal)(IF);
                 (*complexVLocal)(IF) = (*complexV)(IF, iy, ix) + (*complexVLocal)(IF);
                 (*complexWLocal)(IF) = (*complexW)(IF, iy, ix) + (*complexWLocal)(IF);

                 for( int iz = istf; iz <= iedf; ++iz )
                 {
                     (*divergenceOfUStar)(iz) = ckx * (*complexULocal)(iz) + cky * (*complexVLocal)(iz);
                 }
                 CompactDifferenceFirstDerivativeData -> GetDvarDz( istf, iedf, &(*detaDz)(istf), &(*complexWLocal)(istf), &(*dVardZ)(istf), &(*cR)(istf), &(*rC)(istf, -2) );
                 (*divergenceOfUStar)(IF) = (*divergenceOfUStar)(IF) + (*dVardZ)(IF);
                 double dt = GlobalDataBase::GetDoubleParaFromDB("dt");
                 (*deltaPressure)(IF) = (*divergenceOfUStar)(IF) / dt;

                 PoissonSolverData -> PoissonSolverNeumannBC( nSolve, &(*coef_LHS)(istf+1, -2), &(*deltaPressure)(istf), &(*lowerBound)(istf+1), &(*upperBound)(istf+1), nBoundary_Neumann, cZero, cZero, 
                                                              &(*cR)(istf), &(*coefMatrixPressureNeumann)(istf, -2, iy, ix), &(*boundaryMatrixPressureNeumann)(1, 1, iy, ix) );

                 CompactDifferenceFirstDerivativeData -> GetDvarDz( istf, iedf, &(*detaDz)(istf), &(*deltaPressure)(istf), &(*dVardZ)(istf), &(*cR)(istf), &(*rC)(istf, -2) );
                 for( int iz = istf; iz <= iedf; ++iz )
                 {
                     (*complexWLocal)(iz) = (*complexWLocal)(iz) - dt * (*dVardZ)(iz);
                 }

                 (*complexWLocal)(istf) = cZero;
                 (*complexWLocal)(iedf) = cZero;

                 CompactDifferenceFirstDerivativeData -> GetDvarDz( istf, iedf, &(*detaDz)(istf), &(*complexWLocal)(istf), &(*dVardZ)(istf), &(*cR)(istf), &(*rC)(istf, -2) );

                 for( int iz = istf; iz <= iedf; ++iz )
                 {
                     (*omegaZ)(iz) = kyy * (*complexULocal)(iz) - kxx * (*complexVLocal)(iz);
                     (*complexULocal)(iz) = ( cOneI * kxx * (*dVardZ)(iz) + kyy * (*omegaZ)(iz) ) / k2;
                     (*complexVLocal)(iz) = ( cOneI * kyy * (*dVardZ)(iz) - kxx * (*omegaZ)(iz) ) / k2;
                 }
                 (*complexULocal)(istf) = cZero;
                 (*complexULocal)(iedf) = cZero;
                 (*complexVLocal)(istf) = cZero;
                 (*complexVLocal)(iedf) = cZero;

                 (*complexPLocal)(IF) = (*complexPLocal)(IF) + (*deltaPressure)(IF) - 1.0 / (*realPHI)(IF) * (*divergenceOfUStar)(IF);

                 PHComplex dWdZLowerBound = (*dVardZ)(istf);
                 PHComplex dWdZUpperBound = (*dVardZ)(iedf);
                 PHComplex divOfDeltaU1LowerBound = (*complexDivDeltaU1LowerBound)(iy ,ix);
                 PHComplex divOfDeltaU1UpperBound = (*complexDivDeltaU1UpperBound)(iy ,ix);
                 PHComplex divOfDeltaU2LowerBound = (*complexDivDeltaU2LowerBound)(iy ,ix);
                 PHComplex divOfDeltaU2UpperBound = (*complexDivDeltaU2UpperBound)(iy ,ix);

                 PHComplex c1 = (dWdZUpperBound * divOfDeltaU2LowerBound - dWdZLowerBound * divOfDeltaU2UpperBound) 
                              / (divOfDeltaU2UpperBound * divOfDeltaU1LowerBound - divOfDeltaU2LowerBound * divOfDeltaU1UpperBound);
                 PHComplex c2 = (dWdZLowerBound * divOfDeltaU1UpperBound - dWdZUpperBound * divOfDeltaU1LowerBound) 
                              / (divOfDeltaU2UpperBound * divOfDeltaU1LowerBound - divOfDeltaU2LowerBound * divOfDeltaU1UpperBound);

                 for( int iz = istf; iz <= iedf; ++iz )
                 {
                     (*complexULocal)(iz) = (*complexULocal)(iz) + c1 * (*complexDeltaU1)(iz, iy ,ix) + c2 * (*complexDeltaU2)(iz, iy ,ix);
                     (*complexVLocal)(iz) = (*complexVLocal)(iz) + c1 * (*complexDeltaV1)(iz, iy ,ix) + c2 * (*complexDeltaV2)(iz, iy ,ix);
                     (*complexWLocal)(iz) = (*complexWLocal)(iz) + c1 * (*complexDeltaW1)(iz, iy ,ix) + c2 * (*complexDeltaW2)(iz, iy ,ix);
                     (*complexPLocal)(iz) = (*complexPLocal)(iz) + c1 * (*complexDeltaP1)(iz, iy ,ix) + c2 * (*complexDeltaP2)(iz, iy ,ix);
                 }

                 (*complexU)(IF, iy, ix) = (*complexULocal)(IF);
                 (*complexV)(IF, iy, ix) = (*complexVLocal)(IF);
                 (*complexW)(IF, iy, ix) = (*complexWLocal)(IF);
                 (*complexP)(IF, iy, ix) = (*complexPLocal)(IF);                            
            }
        }

        delete realPHI; realPHI = NULL;
        delete complexULocal; complexULocal = NULL;
        delete complexVLocal; complexVLocal = NULL;
        delete complexWLocal; complexWLocal = NULL;
        delete complexPLocal; complexPLocal = NULL;
        delete cR; cR = NULL;
        delete rC; rC = NULL;
        delete dVardZ; dVardZ = NULL;
        delete deltaPressure; deltaPressure = NULL;
        delete divergenceOfUStar; divergenceOfUStar = NULL;
        delete omegaZ; omegaZ = NULL;
    }

    void SpecDiffHybSolver::GetVelocityGradient()
    {
        int ist = GridData -> iStartFourierIndex;
        int ied = GridData -> iEndFourierIndex  ;
        int jst = GridData -> jStartFourierIndex;
        int jed = GridData -> jEndFourierIndex  ;
        int kst = GridData -> kStartFourierIndex;
        int ked = GridData -> kEndFourierIndex;

        Complex1D *varTemp = new Complex1D( Range(ist, ied), fortranArray );
        Complex1D *dvarTemp = new Complex1D( Range(ist, ied), fortranArray );
        Complex1D *workSpace1 = new Complex1D( Range(ist, ied), fortranArray );
        RDouble1D *workSpace2 = new RDouble1D( Range(ist, ied), fortranArray );

        CompactDifferenceFirstDerivativeData -> GetGradientVector( &(*complexU)(ist, jst, kst), &(*complexVelocityGradient)(ist, jst, kst, 1), &(*complexVelocityGradient)(ist, jst, kst, 2), 
                                                                   &(*complexVelocityGradient)(ist, jst, kst, 3), varTemp, dvarTemp, workSpace1, workSpace2, GridData );
        CompactDifferenceFirstDerivativeData -> GetGradientVector( &(*complexV)(ist, jst, kst), &(*complexVelocityGradient)(ist, jst, kst, 4), &(*complexVelocityGradient)(ist, jst, kst, 5), 
                                                                   &(*complexVelocityGradient)(ist, jst, kst, 6), varTemp, dvarTemp, workSpace1, workSpace2, GridData );
        CompactDifferenceFirstDerivativeData -> GetGradientVector( &(*complexW)(ist, jst, kst), &(*complexVelocityGradient)(ist, jst, kst, 7), &(*complexVelocityGradient)(ist, jst, kst, 8), 
                                                                   &(*complexVelocityGradient)(ist, jst, kst, 9), varTemp, dvarTemp, workSpace1, workSpace2, GridData );

        delete varTemp; varTemp = NULL;
        delete dvarTemp; dvarTemp = NULL;
        delete workSpace1; workSpace1 = NULL;
        delete workSpace2; workSpace2 = NULL;
    }

    void SpecDiffHybSolver::GetNonlinearTerm()
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex  ;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex  ;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex;
        Range IF(istf, iedf);
        Range JF(jstf, jedf);
        Range KF(kstf, kedf);

        int istp = GridData -> iStartPhysicalIndex;
        int iedp = GridData -> iEndPhysicalIndex  ;
        int jstp = GridData -> jStartPhysicalIndex;
        int jedp = GridData -> jEndPhysicalIndex  ;
        int kstp = GridData -> kStartPhysicalIndex;
        int kedp = GridData -> kEndPhysicalIndex  ;
        Range IP(istp, iedp);
        Range JP(jstp, jedp);
        Range KP(kstp, kedp);

        int nXYAll = GridData -> nXYAll;
        RDouble rfct_nxny2 = 1.0 / static_cast<double>(nXYAll);

        (*complexVelocityGradient)(IF, JF, KF, 10) = (*complexU)(IF, JF, KF);
        (*complexVelocityGradient)(IF, JF, KF, 11) = (*complexV)(IF, JF, KF);
        (*complexVelocityGradient)(IF, JF, KF, 12) = (*complexW)(IF, JF, KF);
        int nXYZFourier = GridData -> nXYZFourier;
        int nXYZPhysical = GridData -> nXYZPhysical;

        unsigned char opc2r[] = "nff";
        Cp3dfft_btran_c2r_many(reinterpret_cast<double*>(&(*complexVelocityGradient)(istf, jstf, kstf, 1)), nXYZFourier, &(*realVelocityGradient)(istp, jstp, kstp, 1), nXYZPhysical, 12, opc2r);

        (*realVelocitySquare)(IP, JP, KP, 1) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 10); // u*u
        (*realVelocitySquare)(IP, JP, KP, 2) = (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 10); // v*u
        (*realVelocitySquare)(IP, JP, KP, 3) = (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 10); // w*u

        (*realVelocitySquare)(IP, JP, KP, 4) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 11); // u*v
        (*realVelocitySquare)(IP, JP, KP, 5) = (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 11); // v*v
        (*realVelocitySquare)(IP, JP, KP, 6) = (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 11); // w*v

        (*realVelocitySquare)(IP, JP, KP, 7) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 12); // u*w
        (*realVelocitySquare)(IP, JP, KP, 8) = (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 12); // v*w
        (*realVelocitySquare)(IP, JP, KP, 9) = (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 12); // w*w

        (*realVelocitySquare)(IP, JP, KP, 10) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 1)
                                              + (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 2)
                                              + (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 3); // u*du/dx + v*du/dy +w*du/dz
        (*realVelocitySquare)(IP, JP, KP, 11) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 4)
                                              + (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 5)
                                              + (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 6); // u*dv/dx + v*dv/dy +w*dv/dz
        (*realVelocitySquare)(IP, JP, KP, 12) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 7)
                                              + (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 8)
                                              + (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 9); // u*dw/dx + v*dw/dy +w*dw/dz

        unsigned char opr2c[] = "ffn";
        Cp3dfft_ftran_r2c_many(&(*realVelocitySquare)(istp, jstp, kstp, 1), nXYZPhysical, reinterpret_cast<double*>(&(*complexVelocitySquare)(istf, jstf, kstf, 1)), nXYZFourier, 12, opr2c);
        (*complexVelocitySquare)(IF, JF, KF, Range(1, 12)) = (*complexVelocitySquare)(IF, JF, KF, Range(1, 12)) * rfct_nxny2;

        Complex3D *complexDivergenceOfVelocitySquare = new Complex3D(IF, JF, KF, fortranArray);
        Complex1D *varTemp = new Complex1D(IF, fortranArray);
        Complex1D *dVarTemp = new Complex1D(IF, fortranArray);
        Complex3D *dVardZ = new Complex3D(IF, JF, KF, fortranArray); 
        Complex1D *workSpace1 = new Complex1D(IF, fortranArray);
        RDouble1D *workSpace2 = new RDouble1D(IF, fortranArray);
        Complex4D *convectiveTermTemp = new Complex4D(IF, JF, KF, Range(1, 3), fortranArray);

        // ! 0.5*[ vec{u}\dot\grad(u) + \nabla\dot(\vec{u}u) ]
        CompactDifferenceFirstDerivativeData -> GetDivergence( &(*complexVelocitySquare)(istf, jstf, kstf, 1), &(*complexVelocitySquare)(istf, jstf, kstf, 2), &(*complexVelocitySquare)(istf, jstf, kstf, 3),
                                                               complexDivergenceOfVelocitySquare, varTemp, dVarTemp, dVardZ, workSpace1, workSpace2, GridData );
        (*convectiveTermTemp)(IF, JF, KF, 1) = 0.5 * ( (*complexVelocitySquare)(IF, JF, KF, 10) + (*complexDivergenceOfVelocitySquare)(IF, JF, KF) );

        CompactDifferenceFirstDerivativeData -> GetDivergence( &(*complexVelocitySquare)(istf, jstf, kstf, 4), &(*complexVelocitySquare)(istf, jstf, kstf, 5), &(*complexVelocitySquare)(istf, jstf, kstf, 6),
                                                               complexDivergenceOfVelocitySquare, varTemp, dVarTemp, dVardZ, workSpace1, workSpace2, GridData );
        (*convectiveTermTemp)(IF, JF, KF, 2) = 0.5 * ( (*complexVelocitySquare)(IF, JF, KF, 11) + (*complexDivergenceOfVelocitySquare)(IF, JF, KF) );

        CompactDifferenceFirstDerivativeData -> GetDivergence( &(*complexVelocitySquare)(istf, jstf, kstf, 7), &(*complexVelocitySquare)(istf, jstf, kstf, 8), &(*complexVelocitySquare)(istf, jstf, kstf, 9),
                                                               complexDivergenceOfVelocitySquare, varTemp, dVarTemp, dVardZ, workSpace1, workSpace2, GridData );
        (*convectiveTermTemp)(IF, JF, KF, 3) = 0.5 * ( (*complexVelocitySquare)(IF, JF, KF, 12) + (*complexDivergenceOfVelocitySquare)(IF, JF, KF) );
        
        int ivisc = GlobalDataBase::GetIntParaFromDB("ivisc");
        if (ivisc == 0)
        {
            (*realnut)(IF) = 0.0;
            (*realnut_00)(IF) = 0.0;
        }
        else if(ivisc == 1)
        {
            (*realnut)(IF) = 0.0;
            (*realnut_00)(IF) = 0.0;
        }
        else if(ivisc == 2)
        {
            LESData -> Smagorinsky(GridData, realnut, realVelocityGradient);
            PoissonSolverData -> GetCoefMatrixVelocityDirichlet(GridData, DiffBoundaryData, CompactDifferenceSecondDerivativeData, realnut, realnut_00);            
        }
        else if(ivisc == 3)
        {
            LESData -> DynamicSmagorinsky(GridData, CompactDifferenceFirstDerivativeData, realnut, realnut3d, realVelocityGradient, complexU, complexV, complexW);
            PoissonSolverData -> GetCoefMatrixVelocityDirichlet(GridData, DiffBoundaryData, CompactDifferenceSecondDerivativeData, realnut, realnut_00);
        }
        else
        {
            TK_Exit::UnexpectedVarValue( "ivisc", ivisc );
        }

        mode00fix();

        if(ivisc >= 2)
        {
            int inutxyavg = GlobalDataBase::GetIntParaFromDB("inutxyavg");
            if(inutxyavg ==0 || inutxyavg ==1)
            {
                Complex4D *complexDiffTau = LESData -> complexDiffTau;
                (*convectiveTermTemp)(IF, JF, KF, Range(1, 3)) = (*convectiveTermTemp)(IF, JF, KF, Range(1, 3)) + (*complexDiffTau)(IF, JF, KF, Range(1, 3));
            }

            if(inutxyavg == 1 || inutxyavg ==2)
            {
                RDouble1D *detaDz = GridData -> realDetaDz;
                RDouble1D *dnutdz = new RDouble1D(IF, fortranArray);
                RDouble1D *workSpace3 = new RDouble1D(IF, fortranArray);
                CompactDifferenceFirstDerivativeData -> GetDvarDz(istf, iedf, &(*detaDz)(istf), &(*realnut)(istf), &(*dnutdz)(istf), &(*workSpace2)(istf), &(*workSpace3)(istf));

                double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds"); 
                for(int iz = istf; iz <= iedf; ++iz)
                {
                    for(int iy = jstf; iy <= jedf; ++iy)
                    {
                        for(int ix = kstf; ix <= kedf; ++ix)
                        {
                            (*convectiveTermTemp)(iz, iy, ix, 1) = (*convectiveTermTemp)(iz, iy, ix, 1) 
                                - (*dnutdz)(iz) / reynolds * ( (*complexVelocityGradient)(iz, iy, ix, 3) + (*complexVelocityGradient)(iz, iy, ix, 7) );
                            (*convectiveTermTemp)(iz, iy, ix, 2) = (*convectiveTermTemp)(iz, iy, ix, 2) 
                                - (*dnutdz)(iz) / reynolds * ( (*complexVelocityGradient)(iz, iy, ix, 6) + (*complexVelocityGradient)(iz, iy, ix, 8) );
                            (*convectiveTermTemp)(iz, iy, ix, 3) = (*convectiveTermTemp)(iz, iy, ix, 3) 
                                - (*dnutdz)(iz) / reynolds * ( 2.0 * (*complexVelocityGradient)(iz, iy, ix, 9) );
                        }
                    }
                }

                int imode00fix = GlobalDataBase::GetIntParaFromDB("imode00fix");
                using namespace PHMPI;
                int myid = GetCurrentProcessorID();
                if((imode00fix > 0) && (myid ==0))
                {
                    CompactDifferenceFirstDerivativeData -> GetDvarDz(istf, iedf, &(*detaDz)(istf), &(*realnut_00)(istf), &(*dnutdz)(istf), &(*workSpace2)(istf), &(*workSpace3)(istf));

                    int iy = jstf;
                    int ix = kstf;
                    for(int iz = istf; iz <= iedf; ++iz)
                    {
                        (*convectiveTermTemp)(iz, iy, ix, 1) = (*convectiveTermTemp)(iz, iy, ix, 1) 
                            - (*dnutdz)(iz) / reynolds * ( (*complexVelocityGradient)(iz, iy, ix, 3) + (*complexVelocityGradient)(iz, iy, ix, 7) );
                        (*convectiveTermTemp)(iz, iy, ix, 2) = (*convectiveTermTemp)(iz, iy, ix, 2) 
                            - (*dnutdz)(iz) / reynolds * ( (*complexVelocityGradient)(iz, iy, ix, 6) + (*complexVelocityGradient)(iz, iy, ix, 8) );
                        (*convectiveTermTemp)(iz, iy, ix, 3) = (*convectiveTermTemp)(iz, iy, ix, 3) 
                            - (*dnutdz)(iz) / reynolds * ( 2.0 * (*complexVelocityGradient)(iz, iy, ix, 9) );
                    }
                }

                delete dnutdz; dnutdz = NULL;
                delete workSpace3; workSpace3 = NULL;
            }
        }

        (*complexNonlinearTerm)(IF, JF, KF, Range(1, 3)) = (*nonlinearTermCoef)(0) * (*convectiveTermTemp)(IF, JF, KF, Range(1, 3)) + (*nonlinearTermCoef)(-1) * (*convectiveTerm)(IF, JF, KF, Range(1, 3));
        (*convectiveTerm)(IF, JF, KF, Range(1, 3)) = (*convectiveTermTemp)(IF, JF, KF, Range(1, 3));

        delete complexDivergenceOfVelocitySquare; complexDivergenceOfVelocitySquare = NULL;
        delete varTemp; varTemp = NULL;
        delete dVarTemp; dVarTemp = NULL;
        delete dVardZ; dVardZ = NULL;
        delete workSpace1; workSpace1 = NULL;
        delete workSpace2; workSpace2 = NULL;
        delete convectiveTermTemp; convectiveTermTemp = NULL;
    }

    void SpecDiffHybSolver::mode00fix()
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        Range IF(istf, iedf);

        int imode00fix = GlobalDataBase::GetIntParaFromDB("imode00fix");
        if(imode00fix == 0)
        {
            (*realnut_00)(IF) = (*realnut)(IF);
        }
        else
        {
            TK_Exit::UnexpectedVarValue("imode00fix", imode00fix);
        }
    }

    void SpecDiffHybSolver::OutputFieldName()
    {
        fstream file;
        file.open(fieldNameFile.c_str(), ios_base::out);
        if(!file)
        {
            TK_Exit::FileOpenErrorExit(fieldNameFile);
        }
        file << "p\n"
             << "u;velocity\n"
             << "v\n"
             << "w\n";
    }

    void SpecDiffHybSolver::GetViscousTerm()
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex;
        Range IF(istf, iedf);
        Range JF(jstf, jedf);
        Range KF(kstf, kedf);

        Complex1D *c1tmp = new Complex1D(IF, fortranArray);
        Complex1D *c2tmp = new Complex1D(IF, fortranArray);
        Complex3D *c3tmp = new Complex3D(IF, JF, KF, fortranArray);
        Complex3D *cLap_u = new Complex3D(IF, JF, KF, fortranArray);
        Complex3D *cLap_v = new Complex3D(IF, JF, KF, fortranArray);
        Complex3D *cLap_w = new Complex3D(IF, JF, KF, fortranArray);

        CompactDifferenceSecondDerivativeData -> GetLaplacian(GridData, complexU, cLap_u, c1tmp, c2tmp, c3tmp);
        CompactDifferenceSecondDerivativeData -> GetLaplacian(GridData, complexV, cLap_v, c1tmp, c2tmp, c3tmp);
        CompactDifferenceSecondDerivativeData -> GetLaplacian(GridData, complexW, cLap_w, c1tmp, c2tmp, c3tmp);

        for(int iz = istf; iz <= iedf; ++iz)
        {
            (*complexViscousTerm)(iz, JF, KF, 1) = (1.0 + (*realnut)(iz)) * (*cLap_u)(iz, JF, KF);
            (*complexViscousTerm)(iz, JF, KF, 2) = (1.0 + (*realnut)(iz)) * (*cLap_v)(iz, JF, KF);
            (*complexViscousTerm)(iz, JF, KF, 3) = (1.0 + (*realnut)(iz)) * (*cLap_w)(iz, JF, KF);
        }

        int imode00fix = GlobalDataBase::GetIntParaFromDB("imode00fix");
        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        if((imode00fix > 0) && (myid ==0))
        {
            for(int iz = istf; iz <= iedf; ++iz)
            {
                (*complexViscousTerm)(iz, jstf, kstf, 1) = (1.0 + (*realnut)(iz)) * (*cLap_u)(iz, jstf, kstf);
                (*complexViscousTerm)(iz, jstf, kstf, 2) = (1.0 + (*realnut)(iz)) * (*cLap_v)(iz, jstf, kstf);
                (*complexViscousTerm)(iz, jstf, kstf, 3) = (1.0 + (*realnut)(iz)) * (*cLap_w)(iz, jstf, kstf);
            }
        }

        double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");
        (*complexViscousTerm)(IF, JF, KF, Range(1, 3)) = (*complexViscousTerm)(IF, JF, KF, Range(1, 3))/reynolds;

        delete c1tmp; c1tmp = NULL;
        delete c2tmp; c2tmp = NULL;
        delete c3tmp; c3tmp = NULL;
        delete cLap_u; cLap_u = NULL;
        delete cLap_v; cLap_v = NULL;
        delete cLap_w; cLap_w = NULL;
    }

    void SpecDiffHybSolver::VanishingDIVproj()
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex;
        Range IF(istf, iedf);

        RDouble2D *rC = new RDouble2D(IF, Range(-2, 2), fortranArray);
        Complex1D *c1 = new Complex1D(IF, fortranArray);
        Complex1D *c2 = new Complex1D(IF, fortranArray);
        Complex1D *cu_star = new Complex1D(IF, fortranArray);
        Complex1D *cv_star = new Complex1D(IF, fortranArray);
        Complex1D *cw_star = new Complex1D(IF, fortranArray);
        Complex1D *cu_new = new Complex1D(IF, fortranArray);
        Complex1D *cv_new = new Complex1D(IF, fortranArray);
        Complex1D *cw_new = new Complex1D(IF, fortranArray);
        Complex1D *cdudz_star = new Complex1D(IF, fortranArray);
        Complex1D *cdvdz_star = new Complex1D(IF, fortranArray);
        Complex1D *cdwdz = new Complex1D(IF, fortranArray);
        Complex1D *comega_x = new Complex1D(IF, fortranArray);
        Complex1D *comega_y = new Complex1D(IF, fortranArray);
        Complex1D *cR = new Complex1D(IF, fortranArray);

        OutputFieldRestartFile(0);

        const PHComplex cZero = PHComplex(0.0, 0.0);
        const PHComplex cOne = PHComplex(1.0, 0.0);
        int nSolve = iedf - istf;

        RDouble1D *realKX = GridData -> realKX;
        RDouble1D *realKY = GridData -> realKY;
        Complex1D *complexKX = GridData -> complexKX;
        Complex1D *complexKY = GridData -> complexKY;
        RDouble1D *detaDz = GridData -> realDetaDz;

        RDouble2D *coef_LHS = CompactDifferenceSecondDerivativeData -> lHSCoefMatrix;
        RDouble2D *coef_RHS = CompactDifferenceSecondDerivativeData -> rHSCoefMatrix;
        Int1D *lowerBound = CompactDifferenceSecondDerivativeData -> lowerBoundofRHSCoefMatrix;
        Int1D *upperBound = CompactDifferenceSecondDerivativeData -> upperBoundofRHSCoefMatrix;

        Complex2D *complexDivDeltaU1LowerBound = CorrectWallDivergenceData -> complexDivDeltaU1LowerBound;
        Complex2D *complexDivDeltaU1UpperBound = CorrectWallDivergenceData -> complexDivDeltaU1UpperBound;
        Complex2D *complexDivDeltaU2LowerBound = CorrectWallDivergenceData -> complexDivDeltaU2LowerBound;
        Complex2D *complexDivDeltaU2UpperBound = CorrectWallDivergenceData -> complexDivDeltaU2UpperBound;
        Complex3D *complexDeltaU1 = CorrectWallDivergenceData -> complexDeltaU1;
        Complex3D *complexDeltaV1 = CorrectWallDivergenceData -> complexDeltaV1;
        Complex3D *complexDeltaW1 = CorrectWallDivergenceData -> complexDeltaW1;
        Complex3D *complexDeltaU2 = CorrectWallDivergenceData -> complexDeltaU2;
        Complex3D *complexDeltaV2 = CorrectWallDivergenceData -> complexDeltaV2;
        Complex3D *complexDeltaW2 = CorrectWallDivergenceData -> complexDeltaW2;

        for(int ix = kstf; ix <= kedf; ++ix)
        {
            for(int iy = jstf; iy <= jedf; ++iy)
            {
                double rkx = (*realKX)(ix);
                double rky = (*realKY)(iy);
                PHComplex ckx = (*complexKX)(ix);
                PHComplex cky = (*complexKY)(iy);
                double k2  = rkx * rkx + rky * rky;

                if(k2 < TINY)
                {
                    continue;
                }

                (*cu_star)(IF) = (*complexU)(IF, iy, ix);
                (*cv_star)(IF) = (*complexV)(IF, iy, ix);
                (*cw_star)(IF) = (*complexW)(IF, iy, ix);
                CompactDifferenceFirstDerivativeData -> GetDvarDz(istf, iedf, &(*detaDz)(istf), &(*cu_star)(istf), &(*cdudz_star)(istf), &(*cR)(istf), &(*rC)(istf, -2));
                CompactDifferenceFirstDerivativeData -> GetDvarDz(istf, iedf, &(*detaDz)(istf), &(*cv_star)(istf), &(*cdvdz_star)(istf), &(*cR)(istf), &(*rC)(istf, -2));
                for(int iz = istf; iz <= iedf; ++iz)
                {
                    (*comega_x)(iz) = cky * (*cw_star)(iz) - (*cdvdz_star)(iz);
                    (*comega_y)(iz) = -ckx * (*cw_star)(iz) + (*cdudz_star)(iz);
                    (*cw_new)(iz) = cky * (*comega_x)(iz) - ckx * (*comega_y)(iz);
                }
                
                PoissonSolver::PoissonSolverDirichletBC(nSolve, &(*coef_LHS)(istf+1, -2), &(*coef_RHS)(istf+1, -2), &(*cw_new)(istf), k2, &(*lowerBound)(istf+1), &(*upperBound)(istf+1),
                                                        cZero, cZero, &(*rC)(istf, -2), &(*cR)(istf));

                (*cw_new)(istf) = cZero;
                (*cw_new)(iedf) = cZero;
                CorrectWallDivergenceData -> rotProjection();
                
                (*cu_new)(istf) = cZero;
                (*cv_new)(istf) = cZero;
                (*cu_new)(iedf) = cZero;
                (*cv_new)(iedf) = cZero;

                PHComplex cmDIVp = -(*cdwdz)(iedf);
                PHComplex cmDIVm = -(*cdwdz)(istf);
                PHComplex cdd = (*complexDivDeltaU1UpperBound)(iy, ix) * (*complexDivDeltaU2LowerBound)(iy, ix) - (*complexDivDeltaU1LowerBound)(iy, ix) * (*complexDivDeltaU2UpperBound)(iy, ix);
                PHComplex cb1 = ( cmDIVp * (*complexDivDeltaU2LowerBound)(iy, ix) - cmDIVm * (*complexDivDeltaU2UpperBound)(iy, ix)) / cdd;
                PHComplex cb2 = (-cmDIVp * (*complexDivDeltaU1LowerBound)(iy, ix) + cmDIVm * (*complexDivDeltaU1UpperBound)(iy, ix)) / cdd;

                for(int iz = istf; iz <= iedf; ++iz)
                {
                    (*cu_new)(iz) = (*cu_new)(iz) + cb1 * (*complexDeltaU1)(iz, iy, ix) + cb2 * (*complexDeltaU2)(iz, iy, ix);
                    (*cv_new)(iz) = (*cv_new)(iz) + cb1 * (*complexDeltaV1)(iz, iy, ix) + cb2 * (*complexDeltaV2)(iz, iy, ix);
                    (*cw_new)(iz) = (*cw_new)(iz) + cb1 * (*complexDeltaW1)(iz, iy, ix) + cb2 * (*complexDeltaW2)(iz, iy, ix);
                }
            }
        }

        OutputFieldRestartFile(10);

        delete rC; rC = NULL;
        delete c1; c1 = NULL;
        delete c2; c2 = NULL;
        delete cu_star; cu_star = NULL;
        delete cv_star; cv_star = NULL;
        delete cw_star; cw_star = NULL;
        delete cu_new; cu_new = NULL;
        delete cv_new; cv_new = NULL;
        delete cw_new; cw_new = NULL;
        delete cdudz_star; cdudz_star = NULL;
        delete cdvdz_star; cdvdz_star = NULL;
        delete cdwdz; cdwdz = NULL;
        delete comega_x; comega_x = NULL;
        delete comega_y; comega_y = NULL;
        delete cR; cR = NULL;
    }

    void SpecDiffHybSolver::RescaleVelocityField(const int istep)
    {
        int nRescaleFluctuation = GlobalDataBase::GetIntParaFromDB("nRescaleFluctuation");
        int nRescaleUAverage = GlobalDataBase::GetIntParaFromDB("nRescaleUAverage");
        int nFixUAverage = GlobalDataBase::GetIntParaFromDB("nFixUAverage");

        if( (istep < nRescaleFluctuation) || (istep < nRescaleUAverage) )
        {
            RescaleFieldData -> GetUVWMagnitude(GridData, complexU, complexV, complexW);
        }

        if(istep < nRescaleFluctuation)
        {
            RescaleFieldData -> RescaleFluctuation(GridData, complexU, complexV, complexW);
        }

        if(istep < nRescaleUAverage)
        {
            RescaleFieldData -> RescaleUAverage(GridData, complexU);
        }

        if( (nFixUAverage < 0) || (istep < nFixUAverage) )
        {
            RescaleFieldData -> FixUAverage(GridData, complexU);
        }
    }

    void SpecDiffHybSolver::GetStatistics()
    {
        StatisticsData -> GetStatisticsLocal(GridData, complexU, complexV, complexW, complexVelocityGradient);

        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex  ;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex  ;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex;    
        int x0Fall = GridData -> x0Fall;
        int xNFall = GridData -> xNFall;
        int y0Fall = GridData -> y0Fall;
        int yNFall = GridData -> yNFall;

        Range I(istf, iedf);
        Range J(jstf, jedf);
        Range K(kstf, kedf);
        Range XN(x0Fall, xNFall);
        Range YN(y0Fall, yNFall);
        Range N10(1, 10);
        Range N9(1, 9);
        Range N3(1, 3);
        
        RDouble2D *statisticsAllLocal = StatisticsData -> statisticsAllLocal;
        RDouble2D *statisticsAll = StatisticsData -> statisticsAll;
        RDouble2D *epsXYZLocal = StatisticsData -> epsXYZLocal;
        RDouble2D *epsXYZ = StatisticsData -> epsXYZ;

        (*statisticsAll)(I, N10) = (*statisticsAllLocal)(I, N10);
        (*epsXYZ)(I, N9) = (*epsXYZLocal)(I, N9);
        int procIDWithKx0Ky0 = GridData -> procIDWithKx0Ky0;

#ifdef PH_PARALLEL
        int nstat = 10 * (iedf - istf + 1);
        MPI_Reduce(&(*statisticsAllLocal)(istf, 1), &(*statisticsAll)(istf, 1), nstat, MPI_DOUBLE, MPI_SUM, procIDWithKx0Ky0, MPI_COMM_WORLD);

        int neps = 9 * (iedf - istf + 1);
        MPI_Reduce(&(*epsXYZLocal)(istf, 1), &(*epsXYZ)(istf, 1), neps, MPI_DOUBLE, MPI_SUM, procIDWithKx0Ky0, MPI_COMM_WORLD);
#endif

        RDouble3D *energyOfUUKxAll = StatisticsData -> energyOfUUKxAll;
        RDouble3D *energyOfUUKyAll = StatisticsData -> energyOfUUKyAll;
        RDouble3D *energyOfUUKxLocal = StatisticsData -> energyOfUUKxLocal;
        RDouble3D *energyOfUUKyLocal = StatisticsData -> energyOfUUKyLocal;
        Int2D *startAndEndFourierIndexOfAll = GridData -> startAndEndFourierIndexOfAll;
        
        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        int nTProcessor = GetNumberOfProcessor();   

#ifdef PH_PARALLEL
        if(myid == 0)
        {
#endif
            (*energyOfUUKxAll)(I, XN, N3) = 0.0;
            (*energyOfUUKxAll)(I, K, N3) = (*energyOfUUKxLocal)(I, K, N3);

#ifdef PH_PARALLEL
            RDouble3D *energyOfUUKxTmp = new RDouble3D(I, XN, N3, fortranArray);
            for(int ip = 1; ip <= (nTProcessor - 1); ++ip)
            {
                int kstf_ip = (*startAndEndFourierIndexOfAll)(3, ip);
                int kedf_ip = (*startAndEndFourierIndexOfAll)(6, ip);
                int nenergy = (kedf_ip - kstf_ip + 1) * (iedf - istf + 1);
                (*energyOfUUKxTmp)(I, XN, N3) = 0.0;

                for(int in = 1; in <= 3; ++in)
                {
                    MPI_Status status;
                    MPI_Recv(&(*energyOfUUKxTmp)(istf, kstf_ip, in), nenergy, MPI_DOUBLE, ip, ip, MPI_COMM_WORLD, &status);
                    (*energyOfUUKxAll)(I, Range(kstf_ip, kedf_ip), in) = (*energyOfUUKxAll)(I, Range(kstf_ip, kedf_ip), in) + (*energyOfUUKxTmp)(I, Range(kstf_ip, kedf_ip), in);
                }   
            }

            delete energyOfUUKxTmp; energyOfUUKxTmp = NULL;
        }
        else
        {
            int nenergy = (kedf - kstf + 1) * (iedf - istf + 1);
            for(int in = 1; in <= 3; ++in)
            {
                MPI_Send(&(*energyOfUUKxLocal)(istf, kstf, in), nenergy,  MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);
            }
        }
#endif

#ifdef PH_PARALLEL
        if(myid == 0)
        {
#endif
            (*energyOfUUKyAll)(I, YN, N3) = 0.0;
            (*energyOfUUKyAll)(I, J, N3) = (*energyOfUUKyLocal)(I, J, N3);

#ifdef PH_PARALLEL
            RDouble3D *energyOfUUKyTmp = new RDouble3D(I, YN, N3, fortranArray);
            for(int ip = 1; ip <= (nTProcessor - 1); ++ip)
            {
                int jstf_ip = (*startAndEndFourierIndexOfAll)(2, ip);
                int jedf_ip = (*startAndEndFourierIndexOfAll)(5, ip);
                int nenergy = (jedf_ip - jstf_ip + 1) * (iedf - istf + 1);
                (*energyOfUUKyTmp)(I, YN, N3) = 0.0;

                for(int in = 1; in <= 3; ++in)
                {
                    MPI_Status status;
                    MPI_Recv(&(*energyOfUUKyTmp)(istf, jstf_ip, in), nenergy, MPI_DOUBLE, ip, ip, MPI_COMM_WORLD, &status);
                    (*energyOfUUKyAll)(I, Range(jstf_ip, jedf_ip), in) = (*energyOfUUKyAll)(I, Range(jstf_ip, jedf_ip), in) + (*energyOfUUKyTmp)(I, Range(jstf_ip, jedf_ip), in);
                }      
            }

            delete energyOfUUKyTmp; energyOfUUKyTmp = NULL;
        }
        else
        {
            int nenergy = (jedf - jstf + 1) * (iedf - istf + 1);
            for(int in = 1; in <= 3; ++in)
            {
                MPI_Send(&(*energyOfUUKyLocal)(istf, jstf, in), nenergy,  MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);
            }   
        }
#endif

        if(myid != procIDWithKx0Ky0)
        {
            return;
        }

        RDouble1D *uAverageTime = StatisticsData -> uAverageTime;
        RDouble1D *vAverageTime = StatisticsData -> vAverageTime; 
        RDouble1D *wAverageTime = StatisticsData -> wAverageTime;
        RDouble1D *uSquareTime = StatisticsData -> uSquareTime;
        RDouble1D *vSquareTime = StatisticsData -> vSquareTime;
        RDouble1D *wSquareTime = StatisticsData -> wSquareTime; 
        RDouble1D *uMultiplyVTime = StatisticsData -> uMultiplyVTime; 
        RDouble1D *uMultiplyWTime = StatisticsData -> uMultiplyWTime; 
        RDouble1D *vMultiplyWTime = StatisticsData -> vMultiplyWTime;
        RDouble1D *epsilonTime = StatisticsData -> epsilonTime;

        int innz = iedf - istf + 1;
        int iz_half = StatisticsData -> iz_half;
        Range IHalf(istf, iz_half);
        RDouble2D *epsilonXYZTime = new RDouble2D(IHalf, N9, fortranArray);
        RDouble3D *energyOfUUKxT = new RDouble3D(IHalf, XN, N3, fortranArray);
        RDouble3D *energyOfUUKyTT = new RDouble3D(IHalf, YN, N3, fortranArray);
        (*epsilonXYZTime)(IHalf, N9) = 0.0;
        (*energyOfUUKxT)(IHalf, XN, N3) = 0.0;
        (*energyOfUUKyTT)(IHalf, YN, N3) = 0.0;

        for(int iz = istf; iz <= iz_half; ++iz)
        {
            int izup = innz - iz + istf;
            (*uAverageTime)(iz) = 0.5 * ((*statisticsAll)(iz, 1) + (*statisticsAll)(izup, 1));
            (*vAverageTime)(iz) = 0.5 * ((*statisticsAll)(iz, 2) + (*statisticsAll)(izup, 2));
            (*wAverageTime)(iz) = 0.5 * ((*statisticsAll)(iz, 3) + (*statisticsAll)(izup, 3));
            (*uSquareTime)(iz) = 0.5 * ((*statisticsAll)(iz, 4) + (*statisticsAll)(izup, 4));
            (*vSquareTime)(iz) = 0.5 * ((*statisticsAll)(iz, 5) + (*statisticsAll)(izup, 5));
            (*wSquareTime)(iz) = 0.5 * ((*statisticsAll)(iz, 6) + (*statisticsAll)(izup, 6));
            (*uMultiplyVTime)(iz) = 0.5 * ((*statisticsAll)(iz, 7) + (*statisticsAll)(izup, 7));
            (*uMultiplyWTime)(iz) = 0.5 * ((*statisticsAll)(iz, 8) - (*statisticsAll)(izup, 8));
            (*vMultiplyWTime)(iz) = 0.5 * ((*statisticsAll)(iz, 9) - (*statisticsAll)(izup, 9));
            (*epsilonTime)(iz) = 0.5 * ((*statisticsAll)(iz, 10) + (*statisticsAll)(izup, 10));

            (*epsilonXYZTime)(iz, N9) = 0.5 * ((*epsXYZ)(iz, N9) + (*epsXYZ)(izup, N9));

            (*energyOfUUKxT)(iz, XN, N3) = 0.5 * ((*energyOfUUKxAll)(iz, XN, N3) + (*energyOfUUKxAll)(izup, XN, N3));
            (*energyOfUUKyTT)(iz, YN, N3) = 0.5 * ((*energyOfUUKyAll)(iz, YN, N3) + (*energyOfUUKyAll)(izup, YN, N3));
        }

        int nY = GridData -> nY;
        Range YNHalf(y0Fall, nY/2+1);
        RDouble3D *energyOfUUKyT = new RDouble3D(IHalf, YNHalf, N3, fortranArray);
        (*energyOfUUKyT)(IHalf, YNHalf, N3) = 0.0;
        (*energyOfUUKyT)(IHalf, y0Fall, N3) = (*energyOfUUKyTT)(IHalf, y0Fall, N3);
        (*energyOfUUKyT)(IHalf, nY/2+1, N3) = (*energyOfUUKyTT)(IHalf, nY/2+1, N3);
        for(int iy = y0Fall+1; iy <= nY/2; ++iy)
        {
            (*energyOfUUKyT)(IHalf, iy, N3) = (*energyOfUUKyTT)(IHalf, iy, N3) + (*energyOfUUKyTT)(IHalf, nY - iy + 2, N3);
        }

        int istatMeth = GlobalDataBase::GetIntParaFromDB("istatMeth");
        double dt = GlobalDataBase::GetDoubleParaFromDB("dt");
        double tavg = GlobalDataBase::GetDoubleParaFromDB("tavg");
        bool &IsInitStatScratch = StatisticsData -> IsInitStatScratch; 
        RDouble1D *uAverage = StatisticsData -> uAverage;
        RDouble1D *vAverage = StatisticsData -> vAverage;
        RDouble1D *wAverage = StatisticsData -> wAverage;
        RDouble1D *uSquare = StatisticsData -> uSquare; 
        RDouble1D *vSquare = StatisticsData -> vSquare;  
        RDouble1D *wSquare = StatisticsData -> wSquare;  
        RDouble1D *uMultiplyV = StatisticsData -> uMultiplyV;  
        RDouble1D *uMultiplyW = StatisticsData -> uMultiplyW;  
        RDouble1D *vMultiplyW = StatisticsData -> vMultiplyW; 
        RDouble1D *epsilon = StatisticsData -> epsilon;
        RDouble2D *epsilonXYZ = StatisticsData -> epsilonXYZ;
        RDouble3D *energyOfUUKx = StatisticsData -> energyOfUUKx;
        RDouble3D *energyOfUUKy = StatisticsData -> energyOfUUKy;

        if(istatMeth == 0)
        {
            double dtdtavg = dt / tavg;
            if(IsInitStatScratch)
            {
                dtdtavg = 1.0;
            }

            (*uAverage)(IHalf) = (1.0 - dtdtavg) * (*uAverage)(IHalf) + dtdtavg * (*uAverageTime)(IHalf);
            (*vAverage)(IHalf) = (1.0 - dtdtavg) * (*vAverage)(IHalf) + dtdtavg * (*vAverageTime)(IHalf);
            (*wAverage)(IHalf) = (1.0 - dtdtavg) * (*wAverage)(IHalf) + dtdtavg * (*wAverageTime)(IHalf);
            (*uSquare)(IHalf) = (1.0 - dtdtavg) * (*uSquare)(IHalf) + dtdtavg * (*uSquareTime)(IHalf);
            (*vSquare)(IHalf) = (1.0 - dtdtavg) * (*vSquare)(IHalf) + dtdtavg * (*vSquareTime)(IHalf);
            (*wSquare)(IHalf) = (1.0 - dtdtavg) * (*wSquare)(IHalf) + dtdtavg * (*wSquareTime)(IHalf);
            (*uMultiplyV)(IHalf) = (1.0 - dtdtavg) * (*uMultiplyV)(IHalf) + dtdtavg * (*uMultiplyVTime)(IHalf);
            (*uMultiplyW)(IHalf) = (1.0 - dtdtavg) * (*uMultiplyW)(IHalf) + dtdtavg * (*uMultiplyWTime)(IHalf);
            (*vMultiplyW)(IHalf) = (1.0 - dtdtavg) * (*vMultiplyW)(IHalf) + dtdtavg * (*vMultiplyWTime)(IHalf);
            (*epsilon)(IHalf) = (1.0 - dtdtavg) * (*epsilon)(IHalf) + dtdtavg * (*epsilonTime)(IHalf);
            (*epsilonXYZ)(IHalf, N9) = (1.0 - dtdtavg) * (*epsilonXYZ)(IHalf, N9) + dtdtavg * (*epsilonXYZTime)(IHalf, N9);
            (*energyOfUUKx)(IHalf, XN, N3) = (1.0 - dtdtavg) * (*energyOfUUKx)(IHalf, XN, N3) + dtdtavg * (*energyOfUUKxT)(IHalf, XN, N3);
            (*energyOfUUKy)(IHalf, YNHalf, N3) = (1.0 - dtdtavg) * (*energyOfUUKy)(IHalf, YNHalf, N3) + dtdtavg * (*energyOfUUKyT)(IHalf, YNHalf, N3);

            double dUdZLowerBound = real((*complexVelocityGradient)(istf, jstf, kstf, 3));
            double dUdZUpperBound = real((*complexVelocityGradient)(iedf, jstf, kstf, 3));
            double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");
            double &tauWAverage = StatisticsData -> tauWAverage;
            double tauW = 0.5 * ( abs(dUdZLowerBound) +abs(dUdZUpperBound) ) / reynolds;
            tauWAverage = (1.0 - dtdtavg) * tauWAverage + dtdtavg * tauW;

            RDouble1D *nutAverage = StatisticsData -> nutAverage;
            int ivisc = GlobalDataBase::GetIntParaFromDB("ivisc");
            int inutavg = StatisticsData -> inutavg;
            if(ivisc >= 2)
            {
                if(inutavg == 0 || IsInitStatScratch)
                {
                    (*nutAverage)(I) = (*realnut)(I);
                }
                else
                {
                    (*nutAverage)(I) = (1.0 - dtdtavg) * (*nutAverage)(I) + dtdtavg * (*realnut)(I);
                }
                inutavg = 1;
            }

            if(IsInitStatScratch)
            {
                IsInitStatScratch = false;
            }
        }
        else if(istatMeth == 1)
        {           
            int &naddnum = StatisticsData -> naddnum;
            naddnum = naddnum + 1;
            (*uAverage)(IHalf) = (*uAverage)(IHalf) + (*uAverageTime)(IHalf);
            (*vAverage)(IHalf) = (*vAverage)(IHalf) + (*vAverageTime)(IHalf);
            (*wAverage)(IHalf) = (*wAverage)(IHalf) + (*wAverageTime)(IHalf);
            (*uSquare)(IHalf) = (*uSquare)(IHalf) + (*uSquareTime)(IHalf);
            (*vSquare)(IHalf) = (*vSquare)(IHalf) + (*vSquareTime)(IHalf);
            (*wSquare)(IHalf) = (*wSquare)(IHalf) + (*wSquareTime)(IHalf);
            (*uMultiplyV)(IHalf) = (*uMultiplyV)(IHalf) + (*uMultiplyVTime)(IHalf);
            (*uMultiplyW)(IHalf) = (*uMultiplyW)(IHalf) + (*uMultiplyWTime)(IHalf);
            (*vMultiplyW)(IHalf) = (*vMultiplyW)(IHalf) + (*vMultiplyWTime)(IHalf);
            (*epsilon)(IHalf) = (*epsilon)(IHalf) + (*epsilonTime)(IHalf);
            (*epsilonXYZ)(IHalf, N9) = (*epsilonXYZ)(IHalf, N9) + (*epsilonXYZTime)(IHalf, N9);
            (*energyOfUUKx)(IHalf, XN, N3) = (*energyOfUUKx)(IHalf, XN, N3) + (*energyOfUUKxT)(IHalf, XN, N3);
            (*energyOfUUKy)(IHalf, YNHalf, N3) = (*energyOfUUKy)(IHalf, YNHalf, N3) + (*energyOfUUKyT)(IHalf, YNHalf, N3);
        }
        else
        {
            TK_Exit::UnexpectedVarValue("istatMeth", istatMeth);
        }

        delete epsilonXYZTime; epsilonXYZTime = NULL;
        delete energyOfUUKxT; energyOfUUKxT = NULL;
        delete energyOfUUKyTT; energyOfUUKyTT = NULL;
        delete energyOfUUKyT; energyOfUUKyT = NULL;
    }

    void SpecDiffHybSolver::OutputResidual(const int istep)
    {
        string outputdir = "./results";
        outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");
        string residualFile = outputdir + "/monitor.dat";
        string nutFile = outputdir + "/nut.dat";

        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex  ;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex  ;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex; 
        Range I(istf, iedf);
        Range J(jstf, jedf);
        Range K(kstf, kedf);
        (*complexDivergenceU)(I, J, K) = (*complexVelocityGradient)(I, J, K, 1) + (*complexVelocityGradient)(I, J, K, 5) + (*complexVelocityGradient)(I, J, K, 9);

        RDouble1D *rZ = GridData -> realZ;
        double ub = 0.0;
        double tauW = 0.0;
        int procIDWithKx0Ky0 = GridData -> procIDWithKx0Ky0;
        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        if(myid == procIDWithKx0Ky0)
        {
            RDouble1D *uAverage = new RDouble1D(I, fortranArray);
            for(int iz = istf; iz <= iedf; ++iz)
            {
                (*uAverage)(iz) = real((*complexU)(iz, jstf, kstf));
            }
            double dUdZLowerBound = real((*complexVelocityGradient)(istf, jstf, kstf, 3));
            double dUdZUpperBound = real((*complexVelocityGradient)(iedf, jstf, kstf, 3));
            double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");      
            tauW = 0.5 * ( abs(dUdZLowerBound) +abs(dUdZUpperBound) ) / reynolds;

            ub = 0.0;
            for(int iz = (istf+1); iz <= (iedf-1); ++iz)
            {
                ub = ub + 0.5 * ((*rZ)(iz+1) - (*rZ)(iz-1)) * (*uAverage)(iz);
            }
            ub = ub / ((*rZ)(iedf) - (*rZ)(istf));

            delete uAverage; uAverage = NULL;
        }

        RDouble1D *kx = GridData -> realKX;
        RDouble1D *ky = GridData -> realKY;
        double ek = 0.0;
        for(int ix = kstf; ix <= kedf; ++ix)
        {
            for(int iy = jstf; iy <= jedf; ++iy)
            {
                double k2  = (*kx)(ix) * (*kx)(ix) + (*ky)(iy) * (*ky)(iy);
                if( k2 < TINY)
                {
                    continue;
                }
                for(int iz = (istf+1); iz <= (iedf-1); ++iz)
                {
                    ek = ek + 2.0 * real( (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix))
                                        + (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix))
                                        + (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix))
                                        );
                }
            }// !for iy
        }// !for ix
        ek = ek / (iedf - istf - 1);

        double divavg = 0.0;
        for(int iz = istf; iz <= iedf; ++iz)
        {
            double sumDivU = 0.0;
            for(int ix = kstf; ix <= kedf; ++ix)
            {
                for(int iy = jstf; iy <= jedf; ++iy)
                {
                    sumDivU = sumDivU + real( (*complexDivergenceU)(iz, iy, ix) * conj((*complexDivergenceU)(iz, iy, ix)) );
                }
            }
            divavg = divavg + sqrt(2.0 * sumDivU);
        }
        divavg = divavg / (iedf - istf + 1);

#ifdef PH_PARALLEL
        RDouble1D *wmpiloc = new RDouble1D(Range(1, 2), fortranArray);
        RDouble1D *wmpi = new RDouble1D(Range(1, 2), fortranArray);
        (*wmpiloc)(1) = ek;
        (*wmpiloc)(2) = divavg;
        MPI_Reduce(&(*wmpiloc)(1), &(*wmpi)(1), 2, MPI_DOUBLE, MPI_SUM, procIDWithKx0Ky0, MPI_COMM_WORLD);
        ek = (*wmpi)(1);
        divavg = (*wmpi)(2);
        delete wmpiloc; wmpiloc = NULL;
        delete wmpi; wmpi = NULL;
#endif

        if(myid == procIDWithKx0Ky0)
        {
            fstream file;
            file.open(screenOutFile.c_str(), ios_base::out|ios_base::app);
            file << "\n";
            file << "-----------------------------------------------------------------------------------------------------\n";
            file << "istep = " << istep << "\n";
            file << "DIVavg = " << divavg << "\n";
            file.close();
            file.clear();

            double tauWAverage = StatisticsData -> tauWAverage;
            double dt = GlobalDataBase::GetDoubleParaFromDB("dt");           
            file.open(residualFile.c_str(), ios_base::in);
            if ((!file) || IsInitialization)
            {
                file.open(residualFile.c_str(), ios_base::out|ios_base::trunc);
                file << "variables= \n";
                file << "       istep       time        rDIVavg          rUb            rEk            Tauw            rTauw_avg \n";
                file << setw(20) << istep;
                file << setiosflags(ios::right);
                file << setprecision(10);
                file << setiosflags(ios::scientific);
                file << setiosflags(ios::showpoint);
                file << setw(20) << static_cast<double>(istep) * dt 
                     << setw(20) << divavg 
                     << setw(20) << ub 
                     << setw(20) << ek 
                     << setw(20) << tauW 
                     << setw(20) << tauWAverage 
                     << "\n";
                file.close();
                file.clear();
            }
            else
            {
                file.close();
                file.clear();
                file.open(residualFile.c_str(), ios_base::out|ios_base::app);
                file << setw(20) << istep;
                file << setiosflags(ios::right);
                file << setprecision(10);
                file << setiosflags(ios::scientific);
                file << setiosflags(ios::showpoint);
                file << setw(20) << static_cast<double>(istep) * dt 
                     << setw(20) << divavg 
                     << setw(20) << ub 
                     << setw(20) << ek 
                     << setw(20) << tauW 
                     << setw(20) << tauWAverage 
                     << "\n";
                file.close();
                file.clear();
            }
            
            int ivisc = GlobalDataBase::GetIntParaFromDB("ivisc");
            if(ivisc >= 2)
            {
                RDouble1D *nutAverage = StatisticsData -> nutAverage;
                file.open(nutFile.c_str(), ios_base::out);
                if ( !file )
                {
                    TK_Exit::ExceptionExit("could not open nutFile\n");
                }
                file << "variables= z, nut, rnut_avg, rnut_00 \n";
                for(int iz = istf; iz <= iedf; ++iz)
                {
                    file << setiosflags(ios::right);
                    file << setprecision(10);
                    file << setiosflags(ios::scientific);
                    file << setiosflags(ios::showpoint);
                    file << setw(20) << (*rZ)(iz) << setw(20) << (*realnut)(iz) << setw(20) << (*nutAverage)(iz) << setw(20) << (*realnut_00)(iz) << "\n";
                }
            }
        }
    }

    void SpecDiffHybSolver::OutputRestartData(const int istep)
    {
        OutputFieldRestartFile(istep);

        OutputOldTimeData(istep);

        StatisticsData -> OutputStatistics(GridData);
    }

    void SpecDiffHybSolver::OutputOldTimeData(const int istep)
    {
        int nwrest = GlobalDataBase::GetIntParaFromDB("nwrest");
        int inout = istep / nwrest;

        string nonlinearTermFileName, pressureGradientTermFileName, laplacianVelocityTermFileName;
        if(inout % 2  == 0)
        {
            nonlinearTermFileName = nonlinearTermBackupFile;
            pressureGradientTermFileName = pressureGradientTermBackupFile;
            laplacianVelocityTermFileName = laplacianVelocityTermBackupFile;
        }
        else
        {
            nonlinearTermFileName = nonlinearTermFile;
            pressureGradientTermFileName = pressureGradientTermFile;
            laplacianVelocityTermFileName = laplacianVelocityTermFile;
        }

        WriteConvectiveTerm(nonlinearTermFileName);

        int iVortProjection = GlobalDataBase::GetIntParaFromDB("iVortProjection");
        if(iVortProjection != 0)
        {
            TK_Exit::UnexpectedVarValue( "iVortProjection", iVortProjection );
        }

        WriteTimeDependentInformation();
    }

    void SpecDiffHybSolver::WriteConvectiveTerm(const string nonlinearTermFileName)
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex;
        Int1D *fourierIndexSize = GridData -> fourierIndexSize;

        fstream file;
        file.open(nonlinearTermFileName.c_str(), ios_base::out|ios_base::binary);
        int inttmp = 1;
        file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
        file.write(reinterpret_cast<char*>(&(*fourierIndexSize)(1)), 3 * sizeof(int));
        inttmp = 3;
        file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
        int nout = 3 * (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);
        file.write(reinterpret_cast<char*>(&(*convectiveTerm)(istf, jstf, kstf, 1)), nout * sizeof(PHComplex));
        file.close();
        file.clear();
    }

    void SpecDiffHybSolver::WriteTimeDependentInformation()
    {
        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        int procIDWithKx0Ky0 = GridData -> procIDWithKx0Ky0;
        if(myid != procIDWithKx0Ky0)
        {
            return;
        }
        fstream file;
        file.open(timeInformationRestartFile.c_str(), ios_base::out);
        if(!file)
        {
            TK_Exit::FileOpenErrorExit(timeInformationRestartFile);
        }
        file << setw(20) << istep;
        file.close();
        file.clear();
    }

}