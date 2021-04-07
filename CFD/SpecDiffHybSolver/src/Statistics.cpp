#include "Statistics.h"
#include "GlobalDataBase.h"
#include "PHMpi.h"
#include "TK_Exit.h"
#include "PHHeader.h"
#include "TK_Parse.h"

using namespace std;

namespace PHSPACE
{
Statistics::Statistics()
{
}

void Statistics::InitStatisticsData(SpecDiffHybGrid *GridData)
{
    naddnum = 0;
    inutavg = 0;

    int istf = GridData->iStartFourierIndex;
    int iedf = GridData->iEndFourierIndex;
    iz_half = istf + (iedf - istf + 1) / 2;

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData->procIDWithKx0Ky0;

    if (myid == procIDWithKx0Ky0)
    {
        InitStatisticsFileName();
        OutputPlot3dGrid(GridData);
    }

    AllocStatisticsdata(GridData);
    InitStatisticsDataValue(GridData);

    int nstepstat = GlobalDataBase::GetIntParaFromDB("nstepstat");
    if (myid == procIDWithKx0Ky0)
    {
        if (nstepstat < 0)
        {
            InputStatistics(GridData);
        }
        else
        {
            IsInitStatScratch = true;
        }

        WriteStatisticsName();
    }
}

void Statistics::InitStatisticsFileName()
{
    string outputdir = "./results";
    outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");

    zGridPlot3dFile = outputdir + "/Zgrid.x";
    zHalfGridPlot3dFile = outputdir + "/Zgrid_half.x";
    statisticsNameFile = outputdir + "/Statistic.nam";
    statisticsRestartFile = outputdir + "/Statistic.rsta";
    statisticsRestartTimeFile = outputdir + "/Statistic_t.rsta";
    statisticsInformationFile = outputdir + "/Statistic.info";
}

void Statistics::OutputPlot3dGrid(SpecDiffHybGrid *GridData)
{
    int istf = GridData->iStartFourierIndex;
    int iedf = GridData->iEndFourierIndex;
    Int1D *fourierIndexSize = GridData->fourierIndexSize;
    RDouble1D *rZ = GridData->realZ;

    int wordwidth = 20;
    fstream file;
    file.open(zGridPlot3dFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::ExceptionExit("could not open zGridPlot3dFile\n");
    }
    file << setw(wordwidth) << 1 << "\n";
    file << setw(wordwidth) << 1 << setw(wordwidth) << 1 << setw(wordwidth) << (*fourierIndexSize)(1) << "\n";
    for (int iz = istf; iz <= iedf; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << 0.0;
    }
    for (int iz = istf; iz <= iedf; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << 0.0;
    }
    for (int iz = istf; iz <= iedf; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << (*rZ)(iz);
    }
    file.close();
    file.clear();

    int nout = iz_half - istf + 1;

    file.open(zHalfGridPlot3dFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::ExceptionExit("could not open zHalfGridPlot3dFile\n");
    }
    file << setw(wordwidth) << 1 << "\n";
    file << setw(wordwidth) << 1 << setw(wordwidth) << 1 << setw(wordwidth) << nout << "\n";
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << 0.0;
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << 0.0;
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << (*rZ)(iz);
    }
    file.close();
    file.clear();
}

void Statistics::AllocStatisticsdata(SpecDiffHybGrid *GridData)
{
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex  ;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex  ;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex  ;
    int x0Fall = GridData -> x0Fall;
    int xNFall = GridData -> xNFall;
    int y0Fall = GridData -> y0Fall;
    int yNFall = GridData -> yNFall;
    int nY = GridData -> nY;

    Range I(istf, iedf);
    Range IHalf(istf, iz_half);
    Range XN(x0Fall, xNFall);
    Range J(jstf, jedf);
    Range K(kstf, kedf);
    Range YN(y0Fall, yNFall);
    Range YNHalf(y0Fall, nY/2+1);
    Range N10(1, 10);
    Range N9(1, 9);
    Range N3(1, 3);

    uAverage = new RDouble1D(IHalf, fortranArray);
    vAverage = new RDouble1D(IHalf, fortranArray);
    wAverage = new RDouble1D(IHalf, fortranArray);
    uSquare = new RDouble1D(IHalf, fortranArray);
    vSquare = new RDouble1D(IHalf, fortranArray);
    wSquare = new RDouble1D(IHalf, fortranArray);
    uMultiplyV = new RDouble1D(IHalf, fortranArray);
    uMultiplyW = new RDouble1D(IHalf, fortranArray);
    vMultiplyW = new RDouble1D(IHalf, fortranArray);
    epsilon = new RDouble1D(IHalf, fortranArray);

    uAverageTime = new RDouble1D(IHalf, fortranArray);
    vAverageTime = new RDouble1D(IHalf, fortranArray);
    wAverageTime = new RDouble1D(IHalf, fortranArray);
    uSquareTime = new RDouble1D(IHalf, fortranArray);
    vSquareTime = new RDouble1D(IHalf, fortranArray);
    wSquareTime = new RDouble1D(IHalf, fortranArray);
    uMultiplyVTime = new RDouble1D(IHalf, fortranArray);
    uMultiplyWTime = new RDouble1D(IHalf, fortranArray);
    vMultiplyWTime = new RDouble1D(IHalf, fortranArray);
    epsilonTime = new RDouble1D(IHalf, fortranArray);

    statisticsAll = new RDouble2D(I, N10, fortranArray);
    statisticsAllLocal = new RDouble2D(I, N10, fortranArray);
    nutAverage = new RDouble1D(I, fortranArray);

    epsXYZ = new RDouble2D(I, N9, fortranArray); 
    epsXYZLocal = new RDouble2D(I, N9, fortranArray);
    epsilonXYZ = new RDouble2D(IHalf, N9, fortranArray);

    energyOfUUKx = new RDouble3D(IHalf, XN, N3, fortranArray);            
    energyOfUUKxAll = new RDouble3D(I, XN, N3, fortranArray);       
    energyOfUUKxLocal = new RDouble3D(I, K, N3, fortranArray); 
    
    energyOfUUKy = new RDouble3D(IHalf, YNHalf, N3, fortranArray);
    energyOfUUKyAll = new RDouble3D(I, YN, N3, fortranArray);
    energyOfUUKyLocal = new RDouble3D(I, J, N3, fortranArray); 
}

void Statistics::InitStatisticsDataValue(SpecDiffHybGrid *GridData)
{        
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex  ;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex  ;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex  ;
    int x0Fall = GridData -> x0Fall;
    int xNFall = GridData -> xNFall;
    int y0Fall = GridData -> y0Fall;
    int yNFall = GridData -> yNFall;
    int nY = GridData -> nY;

    Range I(istf, iedf);
    Range IHalf(istf, iz_half);
    Range XN(x0Fall, xNFall);
    Range J(jstf, jedf);
    Range K(kstf, kedf);
    Range YN(y0Fall, yNFall);
    Range YNHalf(y0Fall, nY/2+1);
    Range N10(1, 10);
    Range N9(1, 9);
    Range N3(1, 3);

    (*uAverage)(IHalf) = 0.0;
    (*vAverage)(IHalf) = 0.0;
    (*wAverage)(IHalf) = 0.0;
    (*uSquare)(IHalf) = 0.0;
    (*vSquare)(IHalf) = 0.0;
    (*wSquare)(IHalf) = 0.0;
    (*uMultiplyV)(IHalf) = 0.0;
    (*uMultiplyW)(IHalf) = 0.0;
    (*vMultiplyW)(IHalf) = 0.0;
    (*epsilon)(IHalf) = 0.0;

    (*uAverageTime)(IHalf) = 0.0;
    (*vAverageTime)(IHalf) = 0.0;
    (*wAverageTime)(IHalf) = 0.0;
    (*uSquareTime)(IHalf) = 0.0;
    (*vSquareTime)(IHalf) = 0.0;
    (*wSquareTime)(IHalf) = 0.0;
    (*uMultiplyVTime)(IHalf) = 0.0;
    (*uMultiplyWTime)(IHalf) = 0.0;
    (*vMultiplyWTime)(IHalf) = 0.0;
    (*epsilonTime)(IHalf) = 0.0;

    (*statisticsAll)(I, N10) = 0.0;
    (*statisticsAllLocal)(I, N10) = 0.0;
    (*nutAverage)(I) = 0.0;

    (*epsXYZ)(I, N9) = 0.0; 
    (*epsXYZLocal)(I, N9) = 0.0;
    (*epsilonXYZ)(IHalf, N9) = 0.0;

    (*energyOfUUKx)(IHalf, XN, N3) = 0.0;         
    (*energyOfUUKxAll)(I, XN, N3) = 0.0;       
    (*energyOfUUKxLocal)(I, K, N3) = 0.0; 
    
    (*energyOfUUKy)(IHalf, YNHalf, N3) = 0.0;
    (*energyOfUUKyAll)(I, YN, N3) = 0.0;
    (*energyOfUUKyLocal)(I, J, N3) = 0.0;  
}

void Statistics::InputStatistics(SpecDiffHybGrid *GridData)
{
    int instart = GlobalDataBase::GetIntParaFromDB("instart");
    if (instart == 0)
    {
        IsInitStatScratch = true;
        return;
    }

    ReadStatisticsInformation();
    ReadStatisticsRestart(GridData);
}

void Statistics::WriteStatisticsName()
{
    fstream file;
    file.open(statisticsNameFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(statisticsNameFile);
    }

    file << "Uavg; mean velocity\n";
    file << "Vavg\n";
    file << "Wavg\n";
    file << "<uu>\n";
    file << "<vv>\n";
    file << "<ww>\n";
    file << "<uv>\n";
    file << "<uw>\n";
    file << "<vw>\n";
    file << "epsilon\n";
}

void Statistics::ReadStatisticsInformation()
{
    int istatMethTmp = 0;
    int naddnumTmp = 0;

    string line,word;
    string separator  = " =\t\r\n#$,;\"'";
    fstream file;
    file.open(statisticsInformationFile.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(statisticsInformationFile);
    }

    ReadNewLine(file, line);
    line = FindNextWord(line, word, separator);
    from_string<int>(istatMethTmp, word, std::dec);
    line = FindNextWord(line, word, separator);
    from_string<int>(naddnumTmp, word, std::dec);

    int istatMeth = GlobalDataBase::GetIntParaFromDB("istatMeth");
    if (istatMethTmp != istatMeth)
    {
        TK_Exit::ExceptionExit("istatMeth in statisticsInformationFile is NOT equal to that in para.inp!");
    }
    naddnum = naddnumTmp;
    file.close();
    file.clear();
}

void Statistics::ReadStatisticsRestart(SpecDiffHybGrid *GridData)
{
    int istf = GridData -> iStartFourierIndex;

    string line,word;
    string separator = " =\t\r\n#$,;\"'";
    fstream file;
    file.open(statisticsRestartFile.c_str(), ios_base::in);
    if (!file)
    {
        cout << "Warning!!!" << statisticsRestartFile << "does not exist or end during reading!" << "The statisic will begin from scratch!" << "\n" << endl;
        IsInitStatScratch = true;
        file.close();
        file.clear();
    }
    else
    {
        ReadNewLine(file, line);
        ReadNewLine(file, line);

        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*uAverage)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*vAverage)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*wAverage)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*uSquare)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*vSquare)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*wSquare)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*uMultiplyV)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*uMultiplyW)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*vMultiplyW)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*epsilon)(iz);
        }
        file.close();
        file.clear();

        ReadStatisticsSelected(GridData);
    }
}

void Statistics::ReadStatisticsSelected(SpecDiffHybGrid *GridData)
{
    int istf = GridData -> iStartFourierIndex;
    int x0Fall = GridData -> x0Fall;
    int xNFall = GridData -> xNFall;
    int y0Fall = GridData -> y0Fall;
    int nY = GridData -> nY;
    Range IHalf(istf, iz_half);
    Range XN(x0Fall, xNFall);
    Range YNHalf(y0Fall, nY/2+1);
    Range N9(1, 9);

    double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");
    double dreynolds2 = 1.0 / reynolds / reynolds;
    int itmp = 0;
    double rtmp = 0.0;

    string line,word;
    string separator  = " =\t\r\n#$,;\"'";
    fstream file;
    file.open("statistic_selected.dat", ios_base::in);
    if (!file)
    {
        TK_Exit::ExceptionExit("could not open statistic_selected.dat\n");
    }

    ReadNewLine(file, line);
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        ReadNewLine(file, line);

        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);

        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 1), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 2), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 3), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 4), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 5), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 6), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 7), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 8), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 9), word, std::dec );
    }
    (*epsilonXYZ)(IHalf, N9) = (*epsilonXYZ)(IHalf, N9) / dreynolds2;
    file.close();
    file.clear();

    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        // !Read rEuu_kx
        string cheuuf = "./spectra/rEuu_kx_";
        stringstream strmMyid;
        strmMyid << iz;
        cheuuf = cheuuf +  strmMyid.str() + ".dat";
        file.open(cheuuf.c_str(), ios_base::in);
        if (!file)
        {
        cout << "Error in reading in rEuu_kx!!! rEuu_kx=rEuu_ky=0.0 will be used !!!\n" << endl;
        (*energyOfUUKx)(IHalf, XN) = 0.0;
        (*energyOfUUKy)(IHalf, YNHalf) = 0.0;
        }

        ReadNewLine(file, line);
        ReadNewLine(file, line);
        for (int ix = x0Fall; ix <= xNFall; ++ ix)
        {
            ReadNewLine(file, line);

            line = FindNextWord(line, word, separator);

            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKx)(iz, ix, 1), word, std::dec );
            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKx)(iz, ix, 2), word, std::dec );
            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKx)(iz, ix, 3), word, std::dec );
        }

        file.close();
        file.clear();

        //! Read rEuu_ky.
        cheuuf = "./spectra/rEuu_ky_";
        strmMyid;
        strmMyid << iz;
        cheuuf = cheuuf +  strmMyid.str() + ".dat";
        file.open(cheuuf.c_str(), ios_base::in);
        if (!file)
        {
        cout << "Error in reading in rEuu_ky!!! rEuu_kx=rEuu_ky=0.0 will be used !!!\n" << endl;
        (*energyOfUUKx)(IHalf, XN) = 0.0;
        (*energyOfUUKy)(IHalf, YNHalf) = 0.0;
        }

        ReadNewLine(file, line);
        ReadNewLine(file, line);
        for (int iy = y0Fall; iy <= (nY/2+1); ++ iy)
        {
            ReadNewLine(file, line);

            line = FindNextWord(line, word, separator);

            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKy)(iz, iy, 1), word, std::dec );
            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKy)(iz, iy, 2), word, std::dec );
            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKy)(iz, iy, 3), word, std::dec );
        }

        file.close();
        file.clear();
    }
}

void Statistics::GetStatisticsLocal(SpecDiffHybGrid *GridData, Complex3D *complexU, Complex3D *complexV, Complex3D *complexW, Complex4D *complexVelocityGradient)
{
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex  ;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex  ;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex  ;
    int x0Fall = GridData -> x0Fall;
    int y0Fall = GridData -> y0Fall;
    RDouble1D *kx = GridData -> realKX;
    RDouble1D *ky = GridData -> realKY;

    Range I(istf, iedf);
    Range IHalf(istf, iz_half);
    Range J(jstf, jedf);
    Range K(kstf, kedf);
    Range N10(1, 10);
    Range N9(1, 9);
    Range N3(1, 3);

    (*statisticsAllLocal)(I, N10) = 0.0;
    (*epsXYZLocal)(I, N9) = 0.0;
    (*energyOfUUKxLocal)(I, K, N3) = 0.0;
    (*energyOfUUKyLocal)(I, J, N3) = 0.0;

    for (int iz = istf; iz <= iedf; ++ iz)
    {
        if (kstf == x0Fall)
        {
            // !process with (0,:)modes
            int iy0 = jstf;
            if (jstf == y0Fall)
            {
                // !process with (0,0) mode
                iy0 = jstf + 1;
                (*statisticsAllLocal)(iz, 1) = real((*complexU)(iz, jstf, kstf));
            }

            // !Reynolds stress
            PHComplex sumUU = PHComplex(0.0, 0.0);
            PHComplex sumVV = PHComplex(0.0, 0.0);
            PHComplex sumWW = PHComplex(0.0, 0.0);
            PHComplex sumUV = PHComplex(0.0, 0.0);
            PHComplex sumUW = PHComplex(0.0, 0.0);
            PHComplex sumVW = PHComplex(0.0, 0.0);
            for (int iy = iy0; iy <= jedf; ++ iy)
            {
                sumUU = sumUU + (*complexU)(iz, iy, kstf) * conj((*complexU)(iz, iy, kstf));
                sumVV = sumVV + (*complexV)(iz, iy, kstf) * conj((*complexV)(iz, iy, kstf));
                sumWW = sumWW + (*complexW)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf));
                sumUV = sumUV + (*complexU)(iz, iy, kstf) * conj((*complexV)(iz, iy, kstf));
                sumUW = sumUW + (*complexU)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf));
                sumVW = sumVW + (*complexV)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf));
            }
            for (int ix = kstf + 1; ix <= kedf; ++ ix)
            {
                for (int iy = jstf; iy <= jedf; ++ iy)
                {
                    sumUU = sumUU + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix));
                    sumVV = sumVV + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumWW = sumWW + 2.0 * (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                    sumUV = sumUV + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumUW = sumUW + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                    sumVW = sumVW + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                }
            }
            (*statisticsAllLocal)(iz, 4) = real(sumUU);
            (*statisticsAllLocal)(iz, 5) = real(sumVV);
            (*statisticsAllLocal)(iz, 6) = real(sumWW);
            (*statisticsAllLocal)(iz, 7) = real(sumUV);
            (*statisticsAllLocal)(iz, 8) = real(sumUW);
            (*statisticsAllLocal)(iz, 9) = real(sumVW);

            // !dissipation
            int ix = kstf;
            for (int iy = iy0; iy <= jedf; ++ iy)
            {
                double kx2 = (*kx)(ix) * (*kx)(ix);
                double ky2 = (*ky)(iy) * (*ky)(iy);
                double k2 = kx2 + ky2;

                double uu = real( (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix)) );
                double vv = real( (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix)) );
                double ww = real( (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix)) );

                PHComplex complexTmp = k2 * (uu + vv + ww)
                                     + (*complexVelocityGradient)(iz, iy, ix, 3) * conj((*complexVelocityGradient)(iz, iy, ix, 3))
                                     + (*complexVelocityGradient)(iz, iy, ix, 6) * conj((*complexVelocityGradient)(iz, iy, ix, 6))
                                     + (*complexVelocityGradient)(iz, iy, ix, 9) * conj((*complexVelocityGradient)(iz, iy, ix, 9));
                (*statisticsAllLocal)(iz, 10) = (*statisticsAllLocal)(iz, 10) + real(complexTmp);

                (*energyOfUUKxLocal)(iz, ix, 1) = (*energyOfUUKxLocal)(iz, ix, 1) + uu;
                (*epsXYZLocal)(iz, 1) = (*epsXYZLocal)(iz, 1) + kx2 * uu;

                (*energyOfUUKxLocal)(iz, ix, 2) = (*energyOfUUKxLocal)(iz, ix, 2) + vv;
                (*epsXYZLocal)(iz, 2) = (*epsXYZLocal)(iz, 2) + kx2 * vv;

                (*energyOfUUKxLocal)(iz, ix, 3) = (*energyOfUUKxLocal)(iz, ix, 3) + ww;
                (*epsXYZLocal)(iz, 3) = (*epsXYZLocal)(iz, 3) + kx2 * ww;

                (*energyOfUUKyLocal)(iz, iy, 1) = (*energyOfUUKyLocal)(iz, iy, 1) + uu;
                (*epsXYZLocal)(iz, 4) = (*epsXYZLocal)(iz, 4) + ky2 * uu;

                (*energyOfUUKyLocal)(iz, iy, 2) = (*energyOfUUKyLocal)(iz, iy, 2) + vv;
                (*epsXYZLocal)(iz, 5) = (*epsXYZLocal)(iz, 5) + ky2 * vv;

                (*energyOfUUKyLocal)(iz, iy, 3) = (*energyOfUUKyLocal)(iz, iy, 3) + ww;
                (*epsXYZLocal)(iz, 6) = (*epsXYZLocal)(iz, 6) + ky2 * ww;

                complexTmp = (*complexVelocityGradient)(iz, iy, ix, 3) * conj((*complexVelocityGradient)(iz, iy, ix, 3));
                (*epsXYZLocal)(iz, 7) = (*epsXYZLocal)(iz, 7) + real(complexTmp);

                complexTmp = (*complexVelocityGradient)(iz, iy, ix, 6) * conj((*complexVelocityGradient)(iz, iy, ix, 6));
                (*epsXYZLocal)(iz, 8) = (*epsXYZLocal)(iz, 8) + real(complexTmp);

                complexTmp = (*complexVelocityGradient)(iz, iy, ix, 9) * conj((*complexVelocityGradient)(iz, iy, ix, 9));
                (*epsXYZLocal)(iz, 9) = (*epsXYZLocal)(iz, 9) + real(complexTmp);
            }

            for (int ix = kstf+1; ix <= kedf; ++ ix)
            {
                for (int iy = jstf; iy <= jedf; ++ iy)
                {
                    double kx2 = (*kx)(ix) * (*kx)(ix);
                    double ky2 = (*ky)(iy) * (*ky)(iy);
                    double k2 = kx2 + ky2;

                    double uu = real( (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix)) );
                    double vv = real( (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix)) );
                    double ww = real( (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix)) );

                    PHComplex complexTmp = k2 * (uu + vv + ww)
                        + (*complexVelocityGradient)(iz, iy, ix, 3) * conj((*complexVelocityGradient)(iz, iy, ix, 3))
                        + (*complexVelocityGradient)(iz, iy, ix, 6) * conj((*complexVelocityGradient)(iz, iy, ix, 6))
                        + (*complexVelocityGradient)(iz, iy, ix, 9) * conj((*complexVelocityGradient)(iz, iy, ix, 9));
                    (*statisticsAllLocal)(iz, 10) = (*statisticsAllLocal)(iz, 10) + 2.0 * real(complexTmp);

                    (*energyOfUUKxLocal)(iz, ix, 1) = (*energyOfUUKxLocal)(iz, ix, 1) + 2.0 * uu;
                    (*epsXYZLocal)(iz, 1) = (*epsXYZLocal)(iz, 1) + 2.0 * kx2 * uu;

                    (*energyOfUUKxLocal)(iz, ix, 2) = (*energyOfUUKxLocal)(iz, ix, 2) + 2.0 * vv;
                    (*epsXYZLocal)(iz, 2) = (*epsXYZLocal)(iz, 2) + 2.0 * kx2 * vv;

                    (*energyOfUUKxLocal)(iz, ix, 3) = (*energyOfUUKxLocal)(iz, ix, 3) + 2.0 * ww;
                    (*epsXYZLocal)(iz, 3) = (*epsXYZLocal)(iz, 3) + 2.0 * kx2 * ww;

                    (*energyOfUUKyLocal)(iz, iy, 1) = (*energyOfUUKyLocal)(iz, iy, 1) + 2.0 * uu;
                    (*epsXYZLocal)(iz, 4) = (*epsXYZLocal)(iz, 4) + 2.0 * ky2 * uu;

                    (*energyOfUUKyLocal)(iz, iy, 2) = (*energyOfUUKyLocal)(iz, iy, 2) + 2.0 * vv;
                    (*epsXYZLocal)(iz, 5) = (*epsXYZLocal)(iz, 5) + 2.0 * ky2 * vv;

                    (*energyOfUUKyLocal)(iz, iy, 3) = (*energyOfUUKyLocal)(iz, iy, 3) + 2.0 * ww;
                    (*epsXYZLocal)(iz, 6) = (*epsXYZLocal)(iz, 6) + 2.0 * ky2 * ww;

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 3) * conj((*complexVelocityGradient)(iz, iy, ix, 3));
                    (*epsXYZLocal)(iz, 7) = (*epsXYZLocal)(iz, 7) + 2.0 * real(complexTmp);

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 6) * conj((*complexVelocityGradient)(iz, iy, ix, 6));
                    (*epsXYZLocal)(iz, 8) = (*epsXYZLocal)(iz, 8) + 2.0 * real(complexTmp);

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 9) * conj((*complexVelocityGradient)(iz, iy, ix, 9));
                    (*epsXYZLocal)(iz, 9) = (*epsXYZLocal)(iz, 9) + 2.0 * real(complexTmp);
                }// !for iy
            }// !for ix             
        }// !if(kstf == x0Fall)
        else
        {
            // !process with kx!=0
            // !Reynolds stress
            PHComplex sumUU = PHComplex(0.0, 0.0);
            PHComplex sumVV = PHComplex(0.0, 0.0);
            PHComplex sumWW = PHComplex(0.0, 0.0);
            PHComplex sumUV = PHComplex(0.0, 0.0);
            PHComplex sumUW = PHComplex(0.0, 0.0);
            PHComplex sumVW = PHComplex(0.0, 0.0);
            for (int ix = kstf; ix <= kedf; ++ ix)
            {
                for (int iy = jstf; iy <= jedf; ++ iy)
                {
                    sumUU = sumUU + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix));
                    sumVV = sumVV + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumWW = sumWW + 2.0 * (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                    sumUV = sumUV + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumUW = sumUW + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                    sumVW = sumVW + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                }
            }
            (*statisticsAllLocal)(iz, 4) = real(sumUU);
            (*statisticsAllLocal)(iz, 5) = real(sumVV);
            (*statisticsAllLocal)(iz, 6) = real(sumWW);
            (*statisticsAllLocal)(iz, 7) = real(sumUV);
            (*statisticsAllLocal)(iz, 8) = real(sumUW);
            (*statisticsAllLocal)(iz, 9) = real(sumVW);

            // !dissipation
            // !epsilon = -(kx^2+ky^2)*( cu*conjg(cu) +cv*conjg(cv) cw*conjg(cw) ) + cdudz*conjg(cdudz) +cdvdz*conjg(cdvdz) +cdwdz*conjg(cdwdz)
            for (int ix = kstf; ix <= kedf; ++ ix)
            {
                for (int iy = jstf; iy <= jedf; ++ iy)
                {
                    
                    double kx2 = (*kx)(ix) * (*kx)(ix);
                    double ky2 = (*ky)(iy) * (*ky)(iy);
                    double k2 = kx2 + ky2;

                    double uu = real( (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix)) );
                    double vv = real( (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix)) );
                    double ww = real( (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix)) );

                    PHComplex complexTmp = k2 * (uu + vv + ww)
                                         + (*complexVelocityGradient)(iz, iy, ix, 3) * conj((*complexVelocityGradient)(iz, iy, ix, 3))
                                         + (*complexVelocityGradient)(iz, iy, ix, 6) * conj((*complexVelocityGradient)(iz, iy, ix, 6))
                                         + (*complexVelocityGradient)(iz, iy, ix, 9) * conj((*complexVelocityGradient)(iz, iy, ix, 9));
                    (*statisticsAllLocal)(iz, 10) = (*statisticsAllLocal)(iz, 10) + 2.0 * real(complexTmp);

                    (*energyOfUUKxLocal)(iz, ix, 1) = (*energyOfUUKxLocal)(iz, ix, 1) + 2.0 * uu;
                    (*epsXYZLocal)(iz, 1) = (*epsXYZLocal)(iz, 1) + 2.0 * kx2 * uu;

                    (*energyOfUUKxLocal)(iz, ix, 2) = (*energyOfUUKxLocal)(iz, ix, 2) + 2.0 * vv;
                    (*epsXYZLocal)(iz, 2) = (*epsXYZLocal)(iz, 2) + 2.0 * kx2 * vv;

                    (*energyOfUUKxLocal)(iz, ix, 3) = (*energyOfUUKxLocal)(iz, ix, 3) + 2.0 * ww;
                    (*epsXYZLocal)(iz, 3) = (*epsXYZLocal)(iz, 3) + 2.0 * kx2 * ww;

                    (*energyOfUUKyLocal)(iz, iy, 1) = (*energyOfUUKyLocal)(iz, iy, 1) + 2.0 * uu;
                    (*epsXYZLocal)(iz, 4) = (*epsXYZLocal)(iz, 4) + 2.0 * ky2 * uu;

                    (*energyOfUUKyLocal)(iz, iy, 2) = (*energyOfUUKyLocal)(iz, iy, 2) + 2.0 * vv;
                    (*epsXYZLocal)(iz, 5) = (*epsXYZLocal)(iz, 5) + 2.0 * ky2 * vv;

                    (*energyOfUUKyLocal)(iz, iy, 3) = (*energyOfUUKyLocal)(iz, iy, 3) + 2.0 * ww;
                    (*epsXYZLocal)(iz, 6) = (*epsXYZLocal)(iz, 6) + 2.0 * ky2 * ww;

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 3) * conj((*complexVelocityGradient)(iz, iy, ix, 3));
                    (*epsXYZLocal)(iz, 7) = (*epsXYZLocal)(iz, 7) + 2.0 * real(complexTmp);

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 6) * conj((*complexVelocityGradient)(iz, iy, ix, 6));
                    (*epsXYZLocal)(iz, 8) = (*epsXYZLocal)(iz, 8) + 2.0 * real(complexTmp);

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 9) * conj((*complexVelocityGradient)(iz, iy, ix, 9));
                    (*epsXYZLocal)(iz, 9) = (*epsXYZLocal)(iz, 9) + 2.0 * real(complexTmp);
                    
                }//! for iy
            }//! for ix
        }//! else
    }//! for iz
}

void Statistics::OutputStatistics(SpecDiffHybGrid *GridData)
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData -> procIDWithKx0Ky0;

    if (myid != procIDWithKx0Ky0)
    {
        return;
    }
    else
    {
        WriteStatisticsInformation();

        WriteStatisticsRestart(GridData);

        WriteStatisticsSelected(GridData);
    }
}

void Statistics::WriteStatisticsInformation()
{
    int istatMeth = GlobalDataBase::GetIntParaFromDB("istatMeth");
    fstream file;
    file.open(statisticsInformationFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(statisticsInformationFile);
    }
    file << setw(20) << istatMeth << setw(20) << naddnum << "\n";
    file.close();
    file.clear();
}

void Statistics::WriteStatisticsRestart(SpecDiffHybGrid *GridData)
{
    int istf = GridData -> iStartFourierIndex;
    int noutput = iz_half - istf + 1;

    fstream file;
    file.open(statisticsRestartFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(statisticsRestartFile);
    }
    file << setw(20) << 1 << "\n";
    file << setw(20) << 1 << setw(20) << 1 << setw(20) << noutput << setw(20) << 10 << "\n";
    file << setiosflags(ios::right);
    file << setprecision(10);
    file << setiosflags(ios::scientific);
    file << setiosflags(ios::showpoint);
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uAverage)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vAverage)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*wAverage)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uSquare)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vSquare)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*wSquare)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uMultiplyV)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uMultiplyW)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vMultiplyW)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*epsilon)(iz);
    }
    file.close();
    file.clear();

    file.open(statisticsRestartTimeFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(statisticsRestartTimeFile);
    }
    file << setw(20) << 1 << "\n";
    file << setw(20) << 1 << setw(20) << 1 << setw(20) << noutput << setw(20) << 10 << "\n";
    file << setiosflags(ios::right);
    file << setprecision(10);
    file << setiosflags(ios::scientific);
    file << setiosflags(ios::showpoint);
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uAverageTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vAverageTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*wAverageTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uSquareTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vSquareTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*wSquareTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uMultiplyVTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uMultiplyWTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vMultiplyWTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*epsilonTime)(iz);
    }
    file.close();
    file.clear();
}

void Statistics::WriteStatisticsSelected(SpecDiffHybGrid *GridData)
{
    RDouble1D *realZ = GridData -> realZ;
    int istf = GridData -> iStartFourierIndex;
    int x0Fall = GridData -> x0Fall;
    int xNFall = GridData -> xNFall;
    int y0Fall = GridData -> y0Fall;
    int nY = GridData -> nY;
    Range IHalf(istf, iz_half);
    Range XN(x0Fall, xNFall);
    Range YNHalf(y0Fall, nY/2+1);
    Range N9(1, 9);

    double reynolds = GlobalDataBase::GetDoubleParaFromDB("reynolds");
    double dreynolds2 = 1.0 / reynolds / reynolds;

    fstream file;
    file.open("statistic_selected.dat", ios_base::out);
    if (!file)
    {
        TK_Exit::ExceptionExit("could not open statistic_selected.dat\n");
    }
    file << "variables = i, z, zplus, uavg, vavg, wavg, ruu, rvv, rww, epsilon, eps_ux, eps_vx, eps_wx, eps_uy, eps_vy, eps_wy, eps_uz, eps_vz, eps_wz \n";       
    int wordwidth = 20;
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(wordwidth) << iz;
        file << setiosflags(ios::right);
        file << setprecision(10);
        file << setiosflags(ios::scientific);
        file << setiosflags(ios::showpoint);
        file << setw(wordwidth) << (*realZ)(iz)
             << setw(wordwidth) << (1.0 + (*realZ)(iz)) * reynolds
             << setw(wordwidth) << (*uAverage)(iz)
             << setw(wordwidth) << (*vAverage)(iz)
             << setw(wordwidth) << (*wAverage)(iz)
             << setw(wordwidth) << (*uSquare)(iz)
             << setw(wordwidth) << (*vSquare)(iz)
             << setw(wordwidth) << (*wSquare)(iz)
             << setw(wordwidth) << (*epsilon)(iz) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 1) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 2) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 3) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 4) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 5) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 6) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 7) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 8) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 9) * dreynolds2
             << "\n";
    }
    file.close();
    file.clear();

    double realXL = GlobalDataBase::GetDoubleParaFromDB("XL");       
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        //! write rEuu_kx.
        string cheuufx = "./spectra/rEuu_kx_";
        stringstream strmiz;
        strmiz << iz;
        cheuufx = cheuufx + strmiz.str() + ".dat";
        file.open(cheuufx.c_str(), ios_base::out);
        if (!file)
        {
            TK_Exit::ExceptionExit("could not open rEuu_kx.dat\n");
        }
        file << "variables = kx, Euu, Evv, Eww \n";
        file << "zone t = \" yplus = " << setw(wordwidth) << (1.0 + (*realZ)(iz)) * reynolds << "\"" << "\n";
        for (int ix = x0Fall; ix <= xNFall; ++ ix)
        {
            double rk_x = static_cast<double>(ix - x0Fall) / realXL;
            file << setiosflags(ios::right);
            file << setprecision(10);
            file << setiosflags(ios::scientific);
            file << setiosflags(ios::showpoint);
            file << setw(wordwidth) << rk_x
                 << setw(wordwidth) << (*energyOfUUKx)(iz, ix, 1)
                 << setw(wordwidth) << (*energyOfUUKx)(iz, ix, 2)
                 << setw(wordwidth) << (*energyOfUUKx)(iz, ix, 3)
                 << "\n";
        }
        file.close();
        file.clear();
    }

    double realYL = GlobalDataBase::GetDoubleParaFromDB("YL");
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        // !write rEuu_ky
        string cheuufy = "./spectra/rEuu_ky_";
        stringstream strmiz;
        strmiz << iz;
        cheuufy = cheuufy +  strmiz.str() + ".dat";
        file.open(cheuufy.c_str(), ios_base::out);
        if (!file)
        {
            TK_Exit::ExceptionExit("could not open rEuu_ky.dat\n");
        }
        file << "variables = ky, Euu, Evv, Eww \n";
        file << "zone t = \" yplus = " << setw(wordwidth) << (1.0 + (*realZ)(iz)) * reynolds << "\"" << "\n";
        for (int iy = y0Fall; iy <= (nY/2+1); ++ iy)
        {
            double rk_y = static_cast<double>(iy - y0Fall) / realYL;
            file << setiosflags(ios::right);
            file << setprecision(10);
            file << setiosflags(ios::scientific);
            file << setiosflags(ios::showpoint);
            file << setw(wordwidth) << rk_y
                 << setw(wordwidth) << (*energyOfUUKy)(iz, iy, 1)
                 << setw(wordwidth) << (*energyOfUUKy)(iz, iy, 2)
                 << setw(wordwidth) << (*energyOfUUKy)(iz, iy, 3)
                 << "\n";
        }
        file.close();
        file.clear();
    }// !for iz
}

Statistics::~Statistics()
{
    delete uAverage; uAverage = NULL;
    delete vAverage; vAverage = NULL;
    delete wAverage; wAverage = NULL;
    delete uSquare; uSquare = NULL;
    delete vSquare; vSquare = NULL;
    delete wSquare; wSquare = NULL;
    delete uMultiplyV; uMultiplyV = NULL;
    delete uMultiplyW; uMultiplyW = NULL;
    delete vMultiplyW; vMultiplyW = NULL;
    delete epsilon; epsilon = NULL;

    delete uAverageTime; uAverageTime = NULL;
    delete vAverageTime; vAverageTime = NULL;
    delete wAverageTime; wAverageTime = NULL;
    delete uSquareTime; uSquareTime = NULL;
    delete vSquareTime; vSquareTime = NULL;
    delete wSquareTime; wSquareTime = NULL;
    delete uMultiplyVTime; uMultiplyVTime = NULL;
    delete uMultiplyWTime; uMultiplyWTime = NULL;
    delete vMultiplyWTime; vMultiplyWTime = NULL;
    delete epsilonTime; epsilonTime = NULL;

    delete statisticsAll; statisticsAll = NULL;
    delete statisticsAllLocal; statisticsAllLocal = NULL;
    delete nutAverage; nutAverage = NULL;

    delete epsXYZ; epsXYZ = NULL;
    delete epsXYZLocal; epsXYZLocal = NULL;
    delete epsilonXYZ; epsilonXYZ = NULL;

    delete energyOfUUKx; energyOfUUKx = NULL;
    delete energyOfUUKxAll; energyOfUUKxAll = NULL;
    delete energyOfUUKxLocal; energyOfUUKxLocal = NULL;

    delete energyOfUUKy; energyOfUUKy = NULL;
    delete energyOfUUKyAll; energyOfUUKyAll = NULL;
    delete energyOfUUKyLocal; energyOfUUKyLocal = NULL;        
}

}