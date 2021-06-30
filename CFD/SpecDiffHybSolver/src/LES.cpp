#include "LES.h"
#include "CompactDifferenceFirstDerivative.h"
#include "GlobalDataBase.h"
#include "PHMpi.h"
#include "TK_Exit.h"
#include "PHHeader.h"
#include "p3dfft.h"

using namespace std;

namespace PHSPACE
{

	LES::LES()
	{
	}

    void LES::InitLESData(SpecDiffHybGrid *GridData)
    {
        AllocLESData(GridData);

        GetWallDistance(GridData);

        GetDelta(GridData);

        int ivisc = GlobalDataBase::GetIntParaFromDB("ivisc");
        if(ivisc == 3)
        {
            SetTestFilter(GridData);
        }
    }

    void LES::AllocLESData(SpecDiffHybGrid *GridData)
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex  ;
        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex  ;
        int kstf = GridData -> kStartFourierIndex;
		int kedf = GridData -> kEndFourierIndex  ;

		int istp = GridData -> iStartPhysicalIndex;
		int iedp = GridData -> iEndPhysicalIndex  ;
		int jstp = GridData -> jStartPhysicalIndex;
		int jedp = GridData -> jEndPhysicalIndex  ;
		int kstp = GridData -> kStartPhysicalIndex;
        int kedp = GridData -> kEndPhysicalIndex  ;

        Range IF(istf, iedf);
        Range JF(jstf, jedf);
        Range KF(kstf, kedf);

        Range IP(istp, iedp);
        Range JP(jstp, jedp);
        Range KP(kstp, kedp);

        wallDistance = new RDouble1D(IF, fortranArray);
        wallDistancePlus = new RDouble1D(IF, fortranArray);
        vanDriestFunc = new RDouble1D(IF, fortranArray);
        delta = new RDouble1D(IF, fortranArray);
        gkx = new RDouble1D(KF, fortranArray);
        gky = new RDouble1D(JF, fortranArray);
        realTauij = new RDouble4D(IP, JP, KP, Range(1, 6));
        complexTauij = new Complex4D(IF, JF, KF, Range(1, 6));
        complexDiffTau = new Complex4D(IF, JF, KF, Range(1,3));

        (*wallDistance)(IF) = 0.0;
        (*wallDistancePlus)(IF) = 0.0;
        (*vanDriestFunc)(IF) = 0.0;
        (*delta)(IF) = 0.0;
        (*gkx)(KF) = 0.0;
        (*gky)(JF) = 0.0;
        (*realTauij)(IP, JP, KP, Range(1, 6)) = 0.0;
        (*complexTauij)(IF, JF, KF, Range(1, 6)) = PHComplex(0.0, 0.0);
        (*complexDiffTau)(IF, JF, KF, Range(1,3)) = PHComplex(0.0, 0.0);

        int ivisc = GlobalDataBase::GetIntParaFromDB("ivisc");
        if(ivisc == 3)
        {
            dynamicCoef = new RDouble1D(IF, fortranArray);
            complexVarFilter = new Complex4D(IF, JF, KF, Range(1, 21), fortranArray);
            realVarFilter = new RDouble4D(IP, JP, KP, Range(1, 21), fortranArray);
            testFilterGradUij = new Complex4D(IF, JF, KF, Range(1, 9), fortranArray);

            (*dynamicCoef)(IF) = 0.0;
            (*complexVarFilter)(IF, JF, KF, Range(1, 21)) = PHComplex(0.0, 0.0);
            (*realVarFilter)(IP, JP, KP, Range(1, 21)) = 0.0;
            (*testFilterGradUij)(IF, JF, KF, Range(1, 9)) = PHComplex(0.0, 0.0);
        }
    }

    void LES::GetWallDistance(SpecDiffHybGrid *GridData)
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        RDouble1D *realZ = GridData -> realZ;
        
        double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");
        double coefAplus = GlobalDataBase::GetDoubleParaFromDB("coefAplus");
        for(int iz = istf; iz <= iedf; ++iz)
        {
            (*wallDistance)(iz) = abs( (*realZ)(iz) - (*realZ)(istf) );
            (*wallDistance)(iz) = min( (*wallDistance)(iz), abs( (*realZ)(iz) - (*realZ)(iedf) ) );
            (*wallDistancePlus)(iz) = (*wallDistance)(iz) * reynolds;
            (*vanDriestFunc)(iz) = 1.0 - exp( -(*wallDistancePlus)(iz) / coefAplus );
        }
       
        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        if(myid == 0)
        {
            string outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");
            string filename = outputdir + "/rd.dat";
            fstream file;
            file.open(filename.c_str(), ios_base::out);
            if ( !file )
            {
                TK_Exit::FileOpenErrorExit(filename);
            }
            file << "variables=z, rd, rdplus, rVanDriestD \n";
            for(int iz = istf; iz <= iedf; ++iz)
            {
                file << setiosflags(ios::right);
                file << setprecision(10);
                file << setiosflags(ios::scientific);
                file << setiosflags(ios::showpoint);
                file << setw(20) << (*realZ)(iz)
                     << setw(20) << (*wallDistance)(iz)
                     << setw(20) << (*wallDistancePlus)(iz)
                     << setw(20) << (*vanDriestFunc)(iz)
                     << "\n";
            }
            file.close();
            file.clear();
        }
    }

    void LES::GetDelta(SpecDiffHybGrid *GridData)
    {
        int nX = GridData -> nX;
        int nY = GridData -> nY;
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        RDouble1D *realZ = GridData -> realZ;

        double realXL = GlobalDataBase::GetDoubleParaFromDB("XL"); 
        double realYL = GlobalDataBase::GetDoubleParaFromDB("YL");

        const double pi = 2.0 * acos(0.0);
        double lx = 2.0 * pi * realXL;
        double ly = 2.0 * pi * realYL;
        double dx = lx / static_cast<double>(nX);
        double dy = ly / static_cast<double>(nY);

        for(int iz = istf + 1; iz <= iedf - 1; ++iz)
        {
            double dz = 0.5 * abs( (*realZ)(iz + 1) - (*realZ)(iz - 1) );
            (*delta)(iz) = pow(dx * dy * dz, 1.0 / 3.0);
        }

        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        if(myid == 0)
        {
            string outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");
            string filename = outputdir + "/delta.dat";
            fstream file;
            file.open(filename.c_str(), ios_base::out);
            if ( !file )
            {
                TK_Exit::FileOpenErrorExit(filename);
            }
            file << "variables=z, delta \n";
            for(int iz = istf; iz <= iedf; ++iz)
            {
                file << setiosflags(ios::right);
                file << setprecision(10);
                file << setiosflags(ios::scientific);
                file << setiosflags(ios::showpoint);
                file << setw(20) << (*realZ)(iz)
                     << setw(20) << (*delta)(iz)
                     << "\n";
            }
            file.close();
            file.clear();
        }
    }

    void LES::SetTestFilter(SpecDiffHybGrid *GridData)
    {
        const double pi = 2.0 * acos(0.0);
        const double filterWidth = 2.0;

        double XL = GlobalDataBase::GetDoubleParaFromDB("XL");
        double YL = GlobalDataBase::GetDoubleParaFromDB("YL");
        double LX = 2.0 * pi * XL;
        double LY = 2.0 * pi * YL;

        int nX = GridData -> nX;
        int nY = GridData -> nY;
        double dX = LX / static_cast<double>(nX);
        double dY = LY / static_cast<double>(nY);

        int kstf = GridData -> kStartFourierIndex;
        int kedf = GridData -> kEndFourierIndex;
        RDouble1D *realKX = GridData -> realKX;
        for(int ix = kstf; ix <= kedf; ++ix)
        {
            double omegaX = (*realKX)(ix) * dX;
            SetFilterType(filterWidth, omegaX, (*gkx)(ix));
        }

        int jstf = GridData -> jStartFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        RDouble1D *realKY = GridData -> realKY;
        for(int iy = jstf; iy <= jedf; ++iy)
        {
            double omegaY = (*realKY)(iy) * dY;
            SetFilterType(filterWidth, omegaY, (*gky)(iy));
        }
    }

    void LES::SetFilterType(const double filterWidth, double omega, double &gk)
    {
        const double pi = 2.0 * acos(0.0);

        int ifilter = GlobalDataBase::GetIntParaFromDB("ifilter");
        if(ifilter == 1)
        {
            gk = 1.0;
            if(abs(omega) >= pi / filterWidth)
            {
                gk = 0.0;
            }
        }
        else if(ifilter == 2)
        {
            gk = (2.0 + cos(omega)) / 3.0;
        }
        else
        {
            TK_Exit::UnexpectedVarValue("ifilter", ifilter);
        }
    }

    void LES::Smagorinsky(SpecDiffHybGrid *GridData, RDouble1D *realnut, RDouble4D *realVelocityGradient)
    {
        int istf = GridData -> iStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;

        int istp = GridData -> iStartPhysicalIndex;
        int jstp = GridData -> jStartPhysicalIndex;
        int kstp = GridData -> kStartPhysicalIndex;
		int iedp = GridData -> iEndPhysicalIndex;		
		int jedp = GridData -> jEndPhysicalIndex;		
        int kedp = GridData -> kEndPhysicalIndex;

        Range IF(istf, iedf);
        RDouble1D *sijMagnitude = new RDouble1D(IF, fortranArray); 
        RDouble1D *sijMagnitude_loc = new RDouble1D(IF, fortranArray);

        (*sijMagnitude_loc)(IF) = 0.0;
        for(int iz = kstp; iz <= kedp; ++iz)
        {
            for(int iy = jstp; iy <= jedp; ++iy)
            {
                for(int ix = istp; ix <= iedp; ++ix)
                {
                    double s11 = (*realVelocityGradient)(ix, iy, iz, 1);
                    double s22 = (*realVelocityGradient)(ix, iy, iz, 5);
                    double s33 = (*realVelocityGradient)(ix, iy, iz, 9);
                    double s12 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 2) + (*realVelocityGradient)(ix, iy, iz, 4) );
                    double s13 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 3) + (*realVelocityGradient)(ix, iy, iz, 7) );
                    double s23 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 6) + (*realVelocityGradient)(ix, iy, iz, 8) );

                    double sijsij = s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
                    (*sijMagnitude_loc)(iz) = (*sijMagnitude_loc)(iz) + sqrt(2.0 * sijsij);
                }
            }
        }

        int nz = iedf - istf + 1;
        (*sijMagnitude)(IF) = (*sijMagnitude_loc)(IF);
#ifdef PH_PARALLEL       
        MPI_Reduce( &(*sijMagnitude_loc)(istf), &(*sijMagnitude)(istf), nz, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
#endif

        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        if(myid == 0)
        {
            int nXYAll = GridData -> nXYAll;
            (*sijMagnitude)(IF) = (*sijMagnitude)(IF) / static_cast<double>(nXYAll);

            double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");
            double coefSmag = GlobalDataBase::GetDoubleParaFromDB("coefSmag");
            for(int iz = istf; iz <= iedf; ++iz)
            {
                (*realnut)(iz) = reynolds * pow( (coefSmag * (*vanDriestFunc)(iz) * (*delta)(iz)), 2.0 ) * (*sijMagnitude)(iz);
            }
        }

#ifdef PH_PARALLEL
        MPI_Bcast(&(*realnut)(istf), nz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

        delete sijMagnitude; sijMagnitude = NULL;
        delete sijMagnitude_loc; sijMagnitude_loc = NULL;
    }

    void LES::DynamicSmagorinsky(SpecDiffHybGrid *GridData, CompactDifferenceFirstDerivative *CompactDifferenceFirstDerivativeData, 
                                 RDouble1D *realnut, RDouble3D *realnut3d, RDouble4D *realVelocityGradient, Complex3D *complexU, Complex3D *complexV, Complex3D *complexW)
    {
        int istf = GridData -> iStartFourierIndex;
        int jstf = GridData -> jStartFourierIndex;
        int kstf = GridData -> kStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        int kedf = GridData -> kEndFourierIndex;

        int istp = GridData -> iStartPhysicalIndex;
        int jstp = GridData -> jStartPhysicalIndex;
        int kstp = GridData -> kStartPhysicalIndex;
		int iedp = GridData -> iEndPhysicalIndex;		
		int jedp = GridData -> jEndPhysicalIndex;		
        int kedp = GridData -> kEndPhysicalIndex;

        Range IF(istf, iedf);
        Range JF(jstf, jedf);
        Range KF(kstf, kedf);
        RDouble1D *sijMagnitude = new RDouble1D(IF, fortranArray); 
        RDouble1D *sijMagnitude_loc = new RDouble1D(IF, fortranArray);

        //! To get |S| in physical space; 
        (*sijMagnitude_loc)(IF) = 0.0;
        for(int iz = kstp; iz <= kedp; ++iz)
        {
            for(int iy = jstp; iy <= jedp; ++iy)
            {
                for(int ix = istp; ix <= iedp; ++ix)
                {
                    double s11 = (*realVelocityGradient)(ix, iy, iz, 1);
                    double s22 = (*realVelocityGradient)(ix, iy, iz, 5);
                    double s33 = (*realVelocityGradient)(ix, iy, iz, 9);
                    double s12 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 2) + (*realVelocityGradient)(ix, iy, iz, 4) );
                    double s13 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 3) + (*realVelocityGradient)(ix, iy, iz, 7) );
                    double s23 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 6) + (*realVelocityGradient)(ix, iy, iz, 8) );

                    double sijsij = s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
                    (*sijMagnitude_loc)(iz) = (*sijMagnitude_loc)(iz) + sqrt(2.0 * sijsij);
                }
            }
        }

        int nXYAll = GridData -> nXYAll;
        int nz = iedf - istf + 1;
        (*sijMagnitude)(IF) = (*sijMagnitude_loc)(IF);
#ifdef PH_PARALLEL       
        MPI_Reduce( &(*sijMagnitude_loc)(istf), &(*sijMagnitude)(istf), nz, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
#endif
        (*sijMagnitude)(IF) = (*sijMagnitude)(IF) / static_cast<double>(nXYAll);
#ifdef PH_PARALLEL
        MPI_Bcast(&(*sijMagnitude)(istf), nz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

        //! To get uiuj, |S|Sij in physical space; 
        //! realVarFilter(ix, iy, iz, 1 -> 6) = uiuj
        //! realVarFilter(ix, iy, iz, 7 -> 12) = |s|sij
        int inutxyavg = GlobalDataBase::GetIntParaFromDB("inutxyavg");
        for(int iz = kstp; iz <= kedp; ++iz)
        {
            for(int iy = jstp; iy <= jedp; ++iy)
            {
                for(int ix = istp; ix <= iedp; ++ix)
                {
                    double u = (*realVelocityGradient)(ix, iy, iz, 10);
                    double v = (*realVelocityGradient)(ix, iy, iz, 11);
                    double w = (*realVelocityGradient)(ix, iy, iz, 12);

                    (*realVarFilter)(ix, iy, iz, 1) = u * u;
                    (*realVarFilter)(ix, iy, iz, 2) = v * v;
                    (*realVarFilter)(ix, iy, iz, 3) = w * w;
                    (*realVarFilter)(ix, iy, iz, 4) = u * v;
                    (*realVarFilter)(ix, iy, iz, 5) = u * w;
                    (*realVarFilter)(ix, iy, iz, 6) = v * w;

                    double s11 = (*realVelocityGradient)(ix, iy, iz, 1);
                    double s22 = (*realVelocityGradient)(ix, iy, iz, 5);
                    double s33 = (*realVelocityGradient)(ix, iy, iz, 9);
                    double s12 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 2) + (*realVelocityGradient)(ix, iy, iz, 4) );
                    double s13 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 3) + (*realVelocityGradient)(ix, iy, iz, 7) );
                    double s23 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 6) + (*realVelocityGradient)(ix, iy, iz, 8) );
                  
                    if((inutxyavg == 0) || (inutxyavg == 1))
                    {
                        double sijsij = s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
                        sijsij = sqrt(2.0 * sijsij);

                        (*realVarFilter)(ix, iy, iz, 7) = sijsij * s11;
                        (*realVarFilter)(ix, iy, iz, 8) = sijsij * s22;
                        (*realVarFilter)(ix, iy, iz, 9) = sijsij * s33;
                        (*realVarFilter)(ix, iy, iz, 10) = sijsij * s12;
                        (*realVarFilter)(ix, iy, iz, 11) = sijsij * s13;
                        (*realVarFilter)(ix, iy, iz, 12) = sijsij * s23;
                    }
                    else if(inutxyavg == 2)
                    {
                        (*realVarFilter)(ix, iy, iz, 7) = (*sijMagnitude)(iz) * s11;
                        (*realVarFilter)(ix, iy, iz, 8) = (*sijMagnitude)(iz) * s22;
                        (*realVarFilter)(ix, iy, iz, 9) = (*sijMagnitude)(iz) * s33;
                        (*realVarFilter)(ix, iy, iz, 10) = (*sijMagnitude)(iz) * s12;
                        (*realVarFilter)(ix, iy, iz, 11) = (*sijMagnitude)(iz) * s13;
                        (*realVarFilter)(ix, iy, iz, 12) = (*sijMagnitude)(iz) * s23;
                    }
                    else
                    {
                        TK_Exit::UnexpectedVarValue( "inutxyavg", inutxyavg );
                    }
                }
            }
        }

        //! Transfer uiuj, |S|Sij to Spectral space;
        int nXYZPhysical = GridData -> nXYZPhysical;
        int nXYZFourier = GridData -> nXYZFourier;
        unsigned char opr2c[] = "ffn";
        Cp3dfft_ftran_r2c_many(&(*realVarFilter)(istp, jstp, kstp, 1), nXYZPhysical, reinterpret_cast<double*>(&(*complexVarFilter)(istf, jstf, kstf, 1)), nXYZFourier, 12, opr2c);
        (*complexVarFilter)(IF, JF, KF, Range(1, 21)) = (*complexVarFilter)(IF, JF, KF, Range(1, 21)) / static_cast<double>(nXYAll);

       //! To get TF(ui) in Spectral Space;
        for(int ix = kstf; ix <= kedf; ++ix)
        {
            for(int iy = jstf; iy <= jedf; ++iy)
            {
                for(int iz = istf; iz <= iedf; ++iz)
                {
                    (*complexVarFilter)(iz, iy, ix, 13) = (*complexU)(iz, iy, ix);
                    (*complexVarFilter)(iz, iy, ix, 14) = (*complexV)(iz, iy, ix);
                    (*complexVarFilter)(iz, iy, ix, 15) = (*complexW)(iz, iy, ix);
                }
            }
        }

        //! To filter in xy plane by fourier filter;
        FilterXYInSpectralSpace(&(*complexVarFilter)(istf, jstf, kstf, 1), 15, GridData);

        //! To get filtered Sij from filtered ui;
        Complex1D *varTemp = new Complex1D(IF, fortranArray);
        Complex1D *dvarTemp = new Complex1D(IF, fortranArray);
        Complex1D *workSpace1 = new Complex1D(IF, fortranArray);
        RDouble1D *workSpace2 = new RDouble1D(IF, fortranArray);
        CompactDifferenceFirstDerivativeData -> GetGradientVector( &(*complexVarFilter)(istf, jstf, kstf, 13), &(*testFilterGradUij)(istf, jstf, kstf, 1), &(*testFilterGradUij)(istf, jstf, kstf, 2), 
                                                                   &(*testFilterGradUij)(istf, jstf, kstf, 3), varTemp, dvarTemp, workSpace1, workSpace2, GridData );
        CompactDifferenceFirstDerivativeData -> GetGradientVector( &(*complexVarFilter)(istf, jstf, kstf, 14), &(*testFilterGradUij)(istf, jstf, kstf, 4), &(*testFilterGradUij)(istf, jstf, kstf, 5), 
                                                                   &(*testFilterGradUij)(istf, jstf, kstf, 6), varTemp, dvarTemp, workSpace1, workSpace2, GridData );
        CompactDifferenceFirstDerivativeData -> GetGradientVector( &(*complexVarFilter)(istf, jstf, kstf, 15), &(*testFilterGradUij)(istf, jstf, kstf, 7), &(*testFilterGradUij)(istf, jstf, kstf, 8), 
                                                                   &(*testFilterGradUij)(istf, jstf, kstf, 9), varTemp, dvarTemp, workSpace1, workSpace2, GridData );

        for(int ix = kstf; ix <= kedf; ++ix)
        {
            for(int iy = jstf; iy <= jedf; ++iy)
            {
                for(int iz = istf; iz <= iedf; ++iz)
                {
                    (*complexVarFilter)(istf, jstf, kstf, 16) = (*testFilterGradUij)(istf, jstf, kstf, 1);
                    (*complexVarFilter)(istf, jstf, kstf, 17) = (*testFilterGradUij)(istf, jstf, kstf, 5);
                    (*complexVarFilter)(istf, jstf, kstf, 18) = (*testFilterGradUij)(istf, jstf, kstf, 9);
                    (*complexVarFilter)(istf, jstf, kstf, 19) = 0.5 * ((*testFilterGradUij)(istf, jstf, kstf, 2) + (*testFilterGradUij)(istf, jstf, kstf, 4));
                    (*complexVarFilter)(istf, jstf, kstf, 20) = 0.5 * ((*testFilterGradUij)(istf, jstf, kstf, 3) + (*testFilterGradUij)(istf, jstf, kstf, 7));
                    (*complexVarFilter)(istf, jstf, kstf, 21) = 0.5 * ((*testFilterGradUij)(istf, jstf, kstf, 6) + (*testFilterGradUij)(istf, jstf, kstf, 8));
                }
            }
        }

        //! Transfer TF(uiuj), TF(|S|Sij), TF(ui), TF(Sij) back to physical space;
        unsigned char opc2r[] = "nff";
        Cp3dfft_btran_c2r_many(reinterpret_cast<double*>(&(*complexVarFilter)(istf, jstf, kstf, 1)), nXYZFourier, &(*realVarFilter)(istp, jstp, kstp, 1), nXYZPhysical, 21, opc2r);

        //! To get filtered |S| in physical space;
        RDouble1D *testFilterSijMagnitude = new RDouble1D(IF, fortranArray);
        if(inutxyavg == 2)
        {
            (*sijMagnitude_loc)(IF) = 0.0;
            for(int iz = kstp; iz <= kedp; ++iz)
            {
                for(int iy = jstp; iy <= jedp; ++iy)
                {
                    for(int ix = istp; ix <= iedp; ++ix)
                    {
                        double s11 = (*realVarFilter)(ix, iy, iz, 16);
                        double s22 = (*realVarFilter)(ix, iy, iz, 17);
                        double s33 = (*realVarFilter)(ix, iy, iz, 18);
                        double s12 = (*realVarFilter)(ix, iy, iz, 19);
                        double s13 = (*realVarFilter)(ix, iy, iz, 20);
                        double s23 = (*realVarFilter)(ix, iy, iz, 21);

                        double sijsij = s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
                        (*sijMagnitude_loc)(iz) = (*sijMagnitude_loc)(iz) + sqrt(2.0 * sijsij);
                    }
                }
            }

            (*testFilterSijMagnitude)(IF) = (*sijMagnitude_loc)(IF);
#ifdef PH_PARALLEL
            MPI_Reduce(&(*sijMagnitude_loc)(istf), &(*testFilterSijMagnitude)(istf), nz, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
            (*testFilterSijMagnitude)(IF) / static_cast<double>(nXYAll);
#ifdef PH_PARALLEL
            MPI_Bcast(&(*testFilterSijMagnitude)(istf), nz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
        }

        //! To get Lij and Mij;
        double filterWidth = GlobalDataBase::GetDoubleParaFromDB("filterWidth");
        RDouble1D *LijMij_loc = new RDouble1D(IF, fortranArray);
        RDouble1D *MijMij_loc = new RDouble1D(IF, fortranArray);
        (*LijMij_loc)(IF) = 0.0;
        (*MijMij_loc)(IF) = 0.0;
        for(int iz = kstp; iz <= kedp; ++iz)
        {
            for(int ix = istp; ix <= iedp; ++ix)
            {
                for(int iy = jstp; iy <= jedp; ++iy)
                {
                    double testFilterU = (*realVarFilter)(ix, iy, iz, 13);
                    double testFilterV = (*realVarFilter)(ix, iy, iz, 14);
                    double testFilterW = (*realVarFilter)(ix, iy, iz, 15);

                    double L11 = (*realVarFilter)(ix, iy, iz, 1) - testFilterU * testFilterU;
                    double L22 = (*realVarFilter)(ix, iy, iz, 2) - testFilterV * testFilterV;
                    double L33 = (*realVarFilter)(ix, iy, iz, 3) - testFilterW * testFilterW;
                    double L12 = (*realVarFilter)(ix, iy, iz, 4) - testFilterU * testFilterV;
                    double L13 = (*realVarFilter)(ix, iy, iz, 5) - testFilterU * testFilterW;
                    double L23 = (*realVarFilter)(ix, iy, iz, 6) - testFilterV * testFilterW;

                    double testFilterS11 = (*realVarFilter)(ix, iy, iz, 16);
                    double testFilterS22 = (*realVarFilter)(ix, iy, iz, 17);
                    double testFilterS33 = (*realVarFilter)(ix, iy, iz, 18);
                    double testFilterS12 = (*realVarFilter)(ix, iy, iz, 19);
                    double testFilterS13 = (*realVarFilter)(ix, iy, iz, 20);
                    double testFilterS23 = (*realVarFilter)(ix, iy, iz, 21);

                    double M11 = 0.0;
                    double M22 = 0.0;
                    double M33 = 0.0;
                    double M12 = 0.0;
                    double M13 = 0.0;
                    double M23 = 0.0;
                    if(inutxyavg == 0 || inutxyavg == 1)
                    {
                        double testFilterSijSij = testFilterS11 * testFilterS11 + testFilterS22 * testFilterS22 + testFilterS33 * testFilterS33 
                                                + 2.0 * (testFilterS12 * testFilterS12 + testFilterS13 * testFilterS13 + testFilterS23 * testFilterS23);
                        testFilterSijSij = sqrt(2.0 * testFilterSijSij);

                        M11 = filterWidth * filterWidth * testFilterSijSij * testFilterS11 - (*realVarFilter)(ix, iy, iz, 7);
                        M22 = filterWidth * filterWidth * testFilterSijSij * testFilterS22 - (*realVarFilter)(ix, iy, iz, 8);
                        M33 = filterWidth * filterWidth * testFilterSijSij * testFilterS33 - (*realVarFilter)(ix, iy, iz, 9);
                        M12 = filterWidth * filterWidth * testFilterSijSij * testFilterS12 - (*realVarFilter)(ix, iy, iz, 10);
                        M13 = filterWidth * filterWidth * testFilterSijSij * testFilterS13 - (*realVarFilter)(ix, iy, iz, 11);
                        M23 = filterWidth * filterWidth * testFilterSijSij * testFilterS23 - (*realVarFilter)(ix, iy, iz, 12);
                    }
                    else if(inutxyavg == 2)
                    {
                        M11 = filterWidth * filterWidth * (*testFilterSijMagnitude)(iz) * testFilterS11 - (*realVarFilter)(ix, iy, iz, 7);
                        M22 = filterWidth * filterWidth * (*testFilterSijMagnitude)(iz) * testFilterS22 - (*realVarFilter)(ix, iy, iz, 8);
                        M33 = filterWidth * filterWidth * (*testFilterSijMagnitude)(iz) * testFilterS33 - (*realVarFilter)(ix, iy, iz, 9);
                        M12 = filterWidth * filterWidth * (*testFilterSijMagnitude)(iz) * testFilterS12 - (*realVarFilter)(ix, iy, iz, 10);
                        M13 = filterWidth * filterWidth * (*testFilterSijMagnitude)(iz) * testFilterS13 - (*realVarFilter)(ix, iy, iz, 11);
                        M23 = filterWidth * filterWidth * (*testFilterSijMagnitude)(iz) * testFilterS23 - (*realVarFilter)(ix, iy, iz, 12);
                    }
                    else
                    {
                        TK_Exit::UnexpectedVarValue("inutxyavg", inutxyavg);
                    }

                    (*LijMij_loc)(iz) = (*LijMij_loc)(iz) + L11 * M11 + L22 * M22 + L33 * M33 + 2.0 * (L12 * M12 + L13 * M13 + L23 * M23);
                    (*MijMij_loc)(iz) = (*MijMij_loc)(iz) + M11 * M11 + M22 * M22 + M33 * M33 + 2.0 * (M12 * M12 + M13 * M13 + M23 * M23);
                }
            }
        }

        RDouble1D *LijMij = new RDouble1D(IF, fortranArray);
        RDouble1D *MijMij = new RDouble1D(IF, fortranArray);
        (*LijMij)(IF) = (*LijMij_loc)(IF);
        (*MijMij)(IF) = (*MijMij_loc)(IF);
#ifdef PH_PARALLEL
        MPI_Reduce(&(*LijMij_loc)(istf), &(*LijMij)(istf), nz, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&(*MijMij_loc)(istf), &(*MijMij)(istf), nz, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif

        //! To get nut in physical space;
        double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");
        using namespace PHMPI;
        int myid = GetCurrentProcessorID();
        if(myid ==0)
        {
            (*dynamicCoef)(IF) = -0.5 * (*LijMij)(IF) / ((*MijMij)(IF) + TINY);
            for(int iz = istf; iz <= iedf; ++iz)
            {
                (*realnut)(iz) = max(reynolds * (*dynamicCoef)(iz) * (*sijMagnitude)(iz), 0.0);
            }
        }
#ifdef PH_PARALLEL
        MPI_Bcast(&(*dynamicCoef)(istf), nz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(*realnut)(istf), nz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

        if(inutxyavg == 0 || inutxyavg == 1)
        {
            for(int iz = kstp; iz <= kedp; ++iz)
            {
                for(int iy = jstp; iy <= jedp; ++iy)
                {
                    for(int ix = istp; ix <= iedp; ++ix)
                    {
                        double s11 = (*realVelocityGradient)(ix, iy, iz, 1);
                        double s22 = (*realVelocityGradient)(ix, iy, iz, 5);
                        double s33 = (*realVelocityGradient)(ix, iy, iz, 9);
                        double s12 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 2) + (*realVelocityGradient)(ix, iy, iz, 4) );
                        double s13 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 3) + (*realVelocityGradient)(ix, iy, iz, 7) );
                        double s23 = 0.5 * ( (*realVelocityGradient)(ix, iy, iz, 6) + (*realVelocityGradient)(ix, iy, iz, 8) );

                        double sijsij = s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
                        sijsij = sqrt(2.0 * sijsij);

                        double nutTemp = max(reynolds * (*dynamicCoef)(iz) * sijsij, 0.0);
                        if(inutxyavg == 0)
                        {
                            (*realnut3d)(ix, iy, iz) = nutTemp;
                        }
                        else if(inutxyavg == 1)
                        {
                            nutTemp = nutTemp - (*realnut)(iz);
                            (*realnut3d)(ix, iy, iz) = nutTemp;
                        }

                        (*realTauij)(ix, iy, iz, 1) = -2.0 * nutTemp * s11 /reynolds;
                        (*realTauij)(ix, iy, iz, 2) = -2.0 * nutTemp * s22 /reynolds;
                        (*realTauij)(ix, iy, iz, 3) = -2.0 * nutTemp * s33 /reynolds;
                        (*realTauij)(ix, iy, iz, 4) = -2.0 * nutTemp * s12 /reynolds;
                        (*realTauij)(ix, iy, iz, 5) = -2.0 * nutTemp * s13 /reynolds;
                        (*realTauij)(ix, iy, iz, 6) = -2.0 * nutTemp * s23 /reynolds;
                    }
                }
            }

            GetDiffTau(GridData, CompactDifferenceFirstDerivativeData);
        }

        delete sijMagnitude; sijMagnitude = NULL;
        delete sijMagnitude_loc; sijMagnitude_loc = NULL;
        delete varTemp; varTemp = NULL;
        delete dvarTemp; dvarTemp = NULL;
        delete workSpace1; workSpace1 = NULL;
        delete workSpace2; workSpace2 = NULL;
        delete testFilterSijMagnitude; testFilterSijMagnitude = NULL;
        delete LijMij_loc; LijMij_loc = NULL;
        delete MijMij_loc; MijMij_loc = NULL;
        delete LijMij; LijMij = NULL;
        delete MijMij; MijMij = NULL;
    }

    void LES::FilterXYInSpectralSpace(PHComplex *varTemp, const int n, SpecDiffHybGrid * GridData)
    {
        int istf = GridData -> iStartFourierIndex;
        int jstf = GridData -> jStartFourierIndex;
        int kstf = GridData -> kStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        int kedf = GridData -> kEndFourierIndex;

        Range IF(istf, iedf);
        Range JF(jstf, jedf);
        Range KF(kstf, kedf);
        Complex4D var(varTemp, IF, JF, KF, Range(1, n), neverDeleteData, fortranArray );
        for(int in = 1; in <= n; ++in)
        {
            for(int ix = kstf; ix <= kedf; ++ix)
            {
                for(int iy = jstf; iy <= jedf; ++iy)
                {
                    for(int iz = kstf; iz <= kedf; ++iz)
                    {
                        var(iz, iy, ix, in) = var(iz, iy, ix, in) * (*gkx)(ix) * (*gky)(iy);
                    }
                }
            }
        }
    }

    void LES::GetDiffTau(SpecDiffHybGrid *GridData, CompactDifferenceFirstDerivative *CompactDifferenceFirstDerivativeData)
    {
        int istf = GridData -> iStartFourierIndex;
        int jstf = GridData -> jStartFourierIndex;
        int kstf = GridData -> kStartFourierIndex;
        int iedf = GridData -> iEndFourierIndex;
        int jedf = GridData -> jEndFourierIndex;
        int kedf = GridData -> kEndFourierIndex;
        Range IF(istf, iedf);
        Range JF(jstf, jedf);
        Range KF(kstf, kedf);

        int istp = GridData -> iStartPhysicalIndex;
        int jstp = GridData -> jStartPhysicalIndex;
        int kstp = GridData -> kStartPhysicalIndex;

        int nXYZPhysical = GridData -> nXYZPhysical;
        int nXYZFourier = GridData -> nXYZFourier;
        int nXYAll = GridData -> nXYAll;
        unsigned char opr2c[] = "ffn";
        Cp3dfft_ftran_r2c_many(&(*realTauij)(istp, jstp, kstp, 1), nXYZPhysical, reinterpret_cast<double*>(&(*complexTauij)(istf, jstf, kstf, 1)), nXYZFourier, 6, opr2c);
        (*complexTauij)(IF, JF, KF, Range(1, 6)) = (*complexTauij)(IF, JF, KF, Range(1, 6)) / static_cast<double>(nXYAll);

        Complex1D *dTau31Dz = new Complex1D(IF, fortranArray);
        Complex1D *dTau32Dz = new Complex1D(IF, fortranArray);
        Complex1D *dTau33Dz = new Complex1D(IF, fortranArray);
        Complex1D *complexWorkSpace = new Complex1D(IF, fortranArray);
        RDouble1D *realWorkSpace = new RDouble1D(IF, fortranArray);
        RDouble1D *detaDz = GridData -> realDetaDz;
        Complex1D *complexKX = GridData -> complexKX;
        Complex1D *complexKY = GridData -> complexKY;
        for(int ix = kstf; ix <= kedf; ++ix)
        {
            for(int iy = jstf; iy <= jedf; ++iy)
            {
                CompactDifferenceFirstDerivativeData -> GetDvarDz(istf, iedf, &(*detaDz)(istf), &(*complexTauij)(istf, iy, ix, 5), &(*dTau31Dz)(istf), &(*complexWorkSpace)(istf), &(*realWorkSpace)(istf));
                CompactDifferenceFirstDerivativeData -> GetDvarDz(istf, iedf, &(*detaDz)(istf), &(*complexTauij)(istf, iy, ix, 6), &(*dTau32Dz)(istf), &(*complexWorkSpace)(istf), &(*realWorkSpace)(istf));
                CompactDifferenceFirstDerivativeData -> GetDvarDz(istf, iedf, &(*detaDz)(istf), &(*complexTauij)(istf, iy, ix, 3), &(*dTau33Dz)(istf), &(*complexWorkSpace)(istf), &(*realWorkSpace)(istf));

                for(int iz = istf; iz <= iedf; ++iz)
                {
                    (*complexDiffTau)(iz, iy, ix, 1) = (*complexKX)(ix) * (*complexTauij)(iz, iy, ix, 1) + (*complexKY)(iy) * (*complexTauij)(iz, iy, ix, 4) + (*dTau31Dz)(iz);
                    (*complexDiffTau)(iz, iy, ix, 2) = (*complexKX)(ix) * (*complexTauij)(iz, iy, ix, 4) + (*complexKY)(iy) * (*complexTauij)(iz, iy, ix, 2) + (*dTau32Dz)(iz);
                    (*complexDiffTau)(iz, iy, ix, 3) = (*complexKX)(ix) * (*complexTauij)(iz, iy, ix, 5) + (*complexKY)(iy) * (*complexTauij)(iz, iy, ix, 6) + (*dTau33Dz)(iz);
                }
            }
        }

        delete dTau31Dz; dTau31Dz = NULL;
        delete dTau32Dz; dTau32Dz = NULL;
        delete dTau33Dz; dTau33Dz = NULL;
        delete complexWorkSpace; complexWorkSpace = NULL;
        delete realWorkSpace; realWorkSpace = NULL;
    }

	LES::~LES()
	{
        delete wallDistance; wallDistance = NULL;
        delete wallDistancePlus; wallDistancePlus = NULL;
        delete vanDriestFunc; vanDriestFunc = NULL;
        delete delta; delta = NULL;
        delete gkx; gkx = NULL;
        delete gky; gky = NULL;
        delete realTauij; realTauij = NULL;
        delete complexTauij; complexTauij = NULL;
        delete complexDiffTau; complexDiffTau = NULL;

        int ivisc = GlobalDataBase::GetIntParaFromDB("ivisc");
        if(ivisc == 3)
        {
            delete dynamicCoef; dynamicCoef = NULL;
            delete complexVarFilter; complexVarFilter = NULL;
            delete realVarFilter; realVarFilter = NULL;
            delete testFilterGradUij; testFilterGradUij = NULL;
        }
	}

}

