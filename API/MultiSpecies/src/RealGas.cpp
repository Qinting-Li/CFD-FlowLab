#include "RealGas.h"
#include "Constants.h"
#include "PHMatrix.h"
#include "SubsonicBC.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
#include "Glb_Dimension.h"
#include "Geo_SimpleBC.h"
#include "TK_Log.h"
#include "PHIO.h"

namespace PHSPACE
{
namespace GAS_SPACE
{

ThermodynamicFunction::ThermodynamicFunction()
{
    number_of_polynomial_coefficient = 0;
    polynomial_coefficient = 0;
    number_of_temperature_interval = 0;
    enthalpyFitCoefficient = 0;
    polynomial_coefficientT2 = 0;
}

ThermodynamicFunction::~ThermodynamicFunction()
{
    delete []enthalpyFitCoefficient;
    DelPointer2(polynomial_coefficient);
    DelPointer2(polynomial_coefficientT2);
}

void ThermodynamicFunction::Init(int number_of_temperature_interval, int number_of_polynomial_coefficient)
{
    this->number_of_temperature_interval   = number_of_temperature_interval;
    this->number_of_polynomial_coefficient = number_of_polynomial_coefficient;
    polynomial_coefficient = NewPointer2<RDouble>(number_of_temperature_interval, number_of_polynomial_coefficient);
    enthalpyFitCoefficient = new RDouble[number_of_polynomial_coefficient];
    polynomial_coefficientT2 = NewPointer2<RDouble>(number_of_temperature_interval + 1, number_of_polynomial_coefficient);
}

void ThermodynamicFunction::ReadPolynomialCoefficient(fstream &file)
{
    string line,word;
    string separator  = " =\t\r\n#$,;\"'";
    string::size_type npos = string::npos;

    for ( int i_interval = 0; i_interval < number_of_temperature_interval; ++ i_interval )
    {
        getline(file, line);

        RDouble *coef = this->GetPolynomialCoefficient( i_interval );

        for ( int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef )
        {
            line = FindNextWord(line, word, separator);
            from_string< RDouble >(coef[ icoef ], word, std::dec);
        }
    }
}

void ThermodynamicFunction::CreatePolynomialCoefficient(RDouble **curveFitData, RDouble **curveFitDataT2, RDouble *enthalpyFitData)
{
    for (int m = 0; m < number_of_temperature_interval; ++ m)
    {
        RDouble *coef = this->GetPolynomialCoefficient(m);

        for (int n = 0; n < number_of_polynomial_coefficient; ++ n)
        {
            coef[n] = curveFitData[m][n];
        }
    }
    //The curve fit data for the two-temperature model.
    for (int m = 0; m < number_of_temperature_interval + 1; ++ m)
    {
        RDouble *coef = this->GetPolynomialCoefficientT2(m);

        for (int n = 0; n < number_of_polynomial_coefficient; ++ n)
        {
            coef[n] = curveFitDataT2[m][n];
        }
    }

    for (int n = 0; n < number_of_polynomial_coefficient; ++ n)
    {
        enthalpyFitCoefficient[n] = enthalpyFitData[n];
    }
}

void ThermodynamicFunction::ReadPolynomialCoefficient(DataContainer *cdata)
{
    for ( int i_interval = 0; i_interval < number_of_temperature_interval; ++ i_interval )
    {
        RDouble *coef = this->GetPolynomialCoefficient( i_interval );

        for ( int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef )
        {
            cdata->read(&coef[ icoef ], sizeof(RDouble));
        }
    }
    //The curve fit data for the two-temperature model.
    for ( int i_interval = 0; i_interval < number_of_temperature_interval + 1; ++ i_interval )
    {
        RDouble *coef = this->GetPolynomialCoefficientT2( i_interval );

        for ( int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef )
        {
            cdata->read(&coef[ icoef ], sizeof(RDouble));
        }
    }

    for ( int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef )
    {
        cdata->read(&enthalpyFitCoefficient[icoef], sizeof(RDouble));
    }
}

void ThermodynamicFunction::WritePolynomialCoefficient(DataContainer *cdata)
{
    for ( int i_interval = 0; i_interval < number_of_temperature_interval; ++ i_interval )
    {
        RDouble *coef = this->GetPolynomialCoefficient( i_interval );

        for ( int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef )
        {
            cdata->append(&coef[ icoef ], sizeof(RDouble));
        }
    }
    //The curve fit data for the two-temperature model.
    for ( int i_interval = 0; i_interval < number_of_temperature_interval + 1; ++ i_interval )
    {
        RDouble *coef = this->GetPolynomialCoefficientT2( i_interval );

        for ( int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef )
        {
            cdata->append(&coef[ icoef ], sizeof(RDouble));
        }
    }

    for ( int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef )
    {
        cdata->append(&enthalpyFitCoefficient[icoef], sizeof(RDouble));
    }
}

ThermodynamicManager::ThermodynamicManager()
{
    thermoDynamicFunction = 0;
    temperatureRange = 0;
    temperatureStart = 0;
}

ThermodynamicManager::~ThermodynamicManager()
{
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        delete thermoDynamicFunction[ispecies];
    }
    delete [] thermoDynamicFunction;
    delete [] temperatureRange;
    delete [] temperatureStart;
}

void ThermodynamicManager::Init(int numberOfSpecies)
{
    this->numberOfSpecies = numberOfSpecies;
    thermoDynamicFunction = new ThermodynamicFunction * [numberOfSpecies];
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        thermoDynamicFunction[ispecies] = new ThermodynamicFunction();
    }
}

void ThermodynamicManager::Read(fstream &file)
{
    string line,word;
    string separator  = " =\t\r\n#$,;\"'";
    string::size_type npos = string::npos;

    //! Read the number of thermodynamic temperature polynomial coefficients of species.
    SkipLines(file,3);
    getline(file,line);

    //! numberofpolynomialcoefficient is the number of polynomial coefficients(a1-a7,or a1-a6).

    line = FindNextWord(line,word,separator);
    from_string< int >(numberOfTemperatureInterval,word,std::dec);

    line = FindNextWord(line,word,separator);
    from_string< int >(numberofpolynomialcoefficient,word,std::dec);

    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        thermoDynamicFunction[ispecies]->Init(numberOfTemperatureInterval, numberofpolynomialcoefficient);
    }

    temperatureRange = new RDouble[numberOfTemperatureInterval+1];

    //! Read the thermodynamic temperature polynomial coefficients,temperature range and so on of species.
    SkipLines(file,4);
    getline(file, line);

    for ( int m = 0; m < numberOfTemperatureInterval + 1; ++ m )
    {
        line = FindNextWord(line,word,separator);
        from_string< RDouble >(temperatureRange[m],word,std::dec);
    }

    SkipLines(file,2);

    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        SkipLines(file,3);
        thermoDynamicFunction[ispecies]->ReadPolynomialCoefficient(file);
    }
}

void ThermodynamicManager::GenerateData(vector<string>&namelist)
{
    ChemkinFitData myChemkinData;
    ChemkinFitData *currentData = new ChemkinFitData();

    myChemkinData.GenerateFullChemkinData();
    myChemkinData.CreateChemkinFitData(namelist, currentData);

    //! To obtain the basic information of interpolation.
    int nSpecies, nInterval, nCoefNum;
    currentData->GetInterpolationInfo(nSpecies, nInterval, nCoefNum);
    numberOfSpecies = nSpecies;
    numberOfTemperatureInterval = nInterval;
    numberofpolynomialcoefficient = nCoefNum;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        thermoDynamicFunction[ispecies]->Init(numberOfTemperatureInterval, numberofpolynomialcoefficient);
    }

    //! To obtain the temperature range of interpolation.
    RDouble *temperatures = currentData->GetTemperatureRange();
    temperatureRange = new RDouble[numberOfTemperatureInterval + 1];
    for (int m = 0; m < numberOfTemperatureInterval + 1; ++ m)
    {
        temperatureRange[m] = temperatures[m];
    }

    //! To obtain the curve fitting data of interpolation of each species.
    RDouble ***curveFittingData = currentData->GetCurveFittingData();
    RDouble ***curveFittingDataT2 = currentData->GetCurveFittingDataT2();
    RDouble **curveFitData = NewPointer2 <RDouble> (numberOfTemperatureInterval, numberofpolynomialcoefficient);
    RDouble **curveFitDataT2 = NewPointer2 <RDouble> (numberOfTemperatureInterval + 1, numberofpolynomialcoefficient);
    RDouble *tmpFitData = new RDouble[numberofpolynomialcoefficient];
    RDouble **enthalpyFitData = currentData->GetEnthalpyFitData();
    RDouble *benckmarkT2 = currentData->GetStartTemperature();
    temperatureStart = new RDouble [numberOfSpecies];

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        temperatureStart[ispecies] = benckmarkT2[ispecies];

        for (int m = 0; m < numberOfTemperatureInterval; ++ m)
        {
            for (int n = 0; n < numberofpolynomialcoefficient; ++ n)
            {
                curveFitData[m][n] = curveFittingData[ispecies][m][n];
            }
        }

        for (int n = 0; n < numberofpolynomialcoefficient; ++ n)
        {
            tmpFitData[n] = enthalpyFitData[ispecies][n];
        }

        for (int m = 0; m < numberOfTemperatureInterval + 1; ++ m)
        {
            for (int n = 0; n < numberofpolynomialcoefficient; ++ n)
            {
                curveFitDataT2[m][n] = curveFittingDataT2[ispecies][m][n];
            }
        }
        thermoDynamicFunction[ispecies]->CreatePolynomialCoefficient(curveFitData, curveFitDataT2, tmpFitData);
    }

    delete []tmpFitData;
    DelPointer2(curveFitData);
    DelPointer2(curveFitDataT2);
    delete currentData;
}

void ThermodynamicManager::Read(DataContainer *cdata)
{
    cdata->read(&numberOfTemperatureInterval, sizeof(int));
    cdata->read(&numberofpolynomialcoefficient, sizeof(int));

    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        thermoDynamicFunction[ispecies]->Init(numberOfTemperatureInterval, numberofpolynomialcoefficient);
    }

    temperatureStart = new RDouble [numberOfSpecies];
    cdata->read(temperatureStart, numberOfSpecies * sizeof(RDouble));

    temperatureRange = new RDouble[numberOfTemperatureInterval + 1];
    cdata->read(temperatureRange, (numberOfTemperatureInterval + 1) * sizeof(RDouble));

    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        thermoDynamicFunction[ispecies]->ReadPolynomialCoefficient(cdata);
    }
}

void ThermodynamicManager::Write(DataContainer *cdata)
{
    cdata->append(&numberOfTemperatureInterval, sizeof(int));
    cdata->append(&numberofpolynomialcoefficient, sizeof(int));
    cdata->append(temperatureStart, numberOfSpecies * sizeof(RDouble));
    cdata->append(temperatureRange, (numberOfTemperatureInterval + 1) * sizeof(RDouble));
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        thermoDynamicFunction[ispecies]->WritePolynomialCoefficient(cdata);
    }
}

RDouble *ThermodynamicManager::GetPolynomialCoefficient(int iSpecies, int iTemperature)
{
    return thermoDynamicFunction[iSpecies]->GetPolynomialCoefficient(iTemperature);
}

RDouble * ThermodynamicManager::GetPolynomialCoefficientT2(int iSpecies, int iTemperature)
{
    return thermoDynamicFunction[iSpecies]->GetPolynomialCoefficientT2(iTemperature);
}

RDouble* ThermodynamicManager::GetEnthalpyFitCoefficient(int iSpecies)
{
    return thermoDynamicFunction[iSpecies]->GetEnthalpyFitCoefficient();
}

void ThermodynamicManager::GetTemperatureRangeIndex(const RDouble &temperature_dimensional, int &indexOfTemperatureRange)
{
    //! Avoid abnormal situations.
    indexOfTemperatureRange = 0;
    if ( temperatureRange[0] >= temperature_dimensional )
    {
        indexOfTemperatureRange = 0;
    }
    else if ( temperatureRange[numberOfTemperatureInterval] <= temperature_dimensional )
    {
        indexOfTemperatureRange = numberOfTemperatureInterval - 1;
    }
    else
    {
        for ( int i = 0; i < numberOfTemperatureInterval; ++ i )
        {
            if ( ( temperatureRange[i] < temperature_dimensional    ) &&
                 ( temperature_dimensional <= temperatureRange[i+1] )   )
            {
                indexOfTemperatureRange = i;
                return;
            }
        }
    }
    return;
}

RDouble ThermodynamicManager::GetTemperatureMin()
{
    return temperatureRange[0];
}

RDouble ThermodynamicManager::GetTemperatureMax()
{
    return temperatureRange[numberOfTemperatureInterval];
}

int ThermodynamicManager::GetNumberofPolynomialCoefficient()
{
    return this->numberofpolynomialcoefficient;
}

RDouble ThermodynamicManager::GetBenchmarkTemperature(int iSpecies)
{
    return temperatureStart[iSpecies];
}

ReactionRate::ReactionRate()
{
    numberOfReaction = 0;
    afr = 0;
    bfr = 0;
    cfr = 0;

    abr = 0;
    bbr = 0;
    cbr = 0;
}

ReactionRate::~ReactionRate()
{
    delete [] afr;
    delete [] bfr;
    delete [] cfr;

    delete [] abr;
    delete [] bbr;
    delete [] cbr;
}

void ReactionRate::Init(int numberOfReaction)
{
    this->numberOfReaction = numberOfReaction;

    afr = new RDouble [ numberOfReaction ];
    bfr = new RDouble [ numberOfReaction ];
    cfr = new RDouble [ numberOfReaction ];

    abr = new RDouble [ numberOfReaction ];
    bbr = new RDouble [ numberOfReaction ];
    cbr = new RDouble [ numberOfReaction ];
}

void ReactionRate::Read(fstream &file)
{
    string line,word;
    string separator  = " =\t\r\n#$,;\"'";
    string::size_type npos = string::npos;

    if ( numberOfReaction <= 0 ) return;

    //! Read the forward and reverse reaction rate coefficients of each chemical reaction.
    SkipLines(file,3);
    for ( int irection = 0; irection < numberOfReaction; ++ irection )
    {
        getline(file,line);
        //! Read irtmp.
        line = FindNextWord(line, word, separator);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (afr[irection], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (bfr[irection], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (cfr[irection], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (abr[irection], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (bbr[irection], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (cbr[irection], word, std::dec);
    }
}

void ReactionRate::Create(RDouble **reationCoef)
{
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        afr[r] = reationCoef[r][0];
        bfr[r] = reationCoef[r][1];
        cfr[r] = reationCoef[r][2];

        abr[r] = reationCoef[r][3];
        bbr[r] = reationCoef[r][4];
        cbr[r] = reationCoef[r][5];
    }
}

void ReactionRate::Read(DataContainer *cdata)
{
    cdata->read(afr, numberOfReaction * sizeof(RDouble));
    cdata->read(bfr, numberOfReaction * sizeof(RDouble));
    cdata->read(cfr, numberOfReaction * sizeof(RDouble));
                                        
    cdata->read(abr, numberOfReaction * sizeof(RDouble));
    cdata->read(bbr, numberOfReaction * sizeof(RDouble));
    cdata->read(cbr, numberOfReaction * sizeof(RDouble));
}

void ReactionRate::Write(DataContainer *cdata)
{
    cdata->append(afr, numberOfReaction * sizeof(RDouble));
    cdata->append(bfr, numberOfReaction * sizeof(RDouble));
    cdata->append(cfr, numberOfReaction * sizeof(RDouble));

    cdata->append(abr, numberOfReaction * sizeof(RDouble));
    cdata->append(bbr, numberOfReaction * sizeof(RDouble));
    cdata->append(cbr, numberOfReaction * sizeof(RDouble));
}

StoichiometricEquation::StoichiometricEquation()
{
    cvp = 0;
    cvn = 0;
    cpz = 0;
    numberOfReaction = 0;
    numberOfSpecies = 0;

    reactionType = 0;
    isCollisionType = 0;
}

StoichiometricEquation::~StoichiometricEquation()
{
    DelPointer2(cvp);
    DelPointer2(cvn);
    DelPointer2(cpz);

    delete []reactionType;
    delete []isCollisionType;
}

void StoichiometricEquation::Init(int numberOfReaction, int numberOfSpecies)
{
    this->numberOfReaction = numberOfReaction;
    this->numberOfSpecies  = numberOfSpecies ;

    reactionType = new int[numberOfReaction];
    isCollisionType = new bool[numberOfReaction];

    cvp = NewPointer2<int>    (numberOfReaction, numberOfSpecies);
    cvn = NewPointer2<int>    (numberOfReaction, numberOfSpecies);
    cpz = NewPointer2<RDouble>    (numberOfReaction, numberOfSpecies);
}

void StoichiometricEquation::Read(fstream &file)
{
    string line,word;
    string separator  = " =\t\r\n#$,;\"'";
    string::size_type npos = string::npos;

    if ( numberOfReaction <= 0 ) return;

    SkipLines(file,4);

    //! Read chemical reaction equation coefficients of species.
    for ( int irection = 0; irection < numberOfReaction; ++ irection )
    {
        getline(file,line);
        //! Read irtmp.
        line = FindNextWord(line,word,separator);

        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            line = FindNextWord(line,word,separator);
            from_string< int >(cvp[irection][ispecies],word,std::dec);
        }
    }

    SkipLines(file,4);
    //! Read chemical reaction equation coefficients of species.
    for ( int irection = 0; irection < numberOfReaction; ++ irection )
    {
        getline(file,line);
        //! Read irtmp.
        line = FindNextWord(line,word,separator);

        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            line = FindNextWord(line,word,separator);
            from_string< int >(cvn[irection][ispecies],word,std::dec);
        }
    }

    SkipLines(file,4);
    //! Read chemical reaction equation coefficients of each collision body.
    for ( int irection = 0; irection < numberOfReaction; ++ irection )
    {
        getline(file,line);
        //! Read irtmp.
        line = FindNextWord(line,word,separator);

        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            line = FindNextWord(line,word,separator);
            from_string< RDouble >(cpz[irection][ispecies], word, std::dec);
        }
    }

    //Read the types of chemical reactions.
    SkipLines(file, 3);
    getline(file,line);
    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        line = FindNextWord(line, word, separator);
        from_string<int>(reactionType[irection], word, std::dec);

        isCollisionType[irection] = 0;
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            //Determine whether the third collision body exists.
            if (cpz[irection][ispecies] > 0.0) 
            {
                isCollisionType[irection] = 1;
                break;
            }
        }
    }
}

void StoichiometricEquation::Create(int **forwardCoef, int **backwardCoef, RDouble **collisionCoef, int **reactionFlag)
{
    //! To create the coefficient of forward reaction equation.
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            cvp[r][s] = forwardCoef[r][s];
        }
    }

    //! To create the coefficient of backward reaction equation.
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            cvn[r][s] = backwardCoef[r][s];
        }
    }

    //! To create the coefficient of the third collision body in the reaction equation.
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            cpz[r][s] = collisionCoef[r][s];
        }
    }

    //! To generate the types of chemical reactions.
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        reactionType[r] = reactionFlag[r][1];
        if (reactionFlag[r][0] == 1)
            isCollisionType[r] = true;
        else
            isCollisionType[r] = false;
    }
}

void StoichiometricEquation::Read(DataContainer *cdata)
{
    for ( int irection = 0; irection < numberOfReaction; ++ irection )
    {
        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            cdata->read(&cvp[irection][ispecies], sizeof(int));
        }
    }

    for ( int irection = 0; irection < numberOfReaction; ++ irection )
    {
        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            cdata->read(&cvn[irection][ispecies], sizeof(int));
        }
    }

    for ( int irection = 0; irection < numberOfReaction; ++ irection )
    {
        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            cdata->read(&cpz[irection][ispecies], sizeof(RDouble));
        }
    }

    cdata->read(reactionType, numberOfReaction * sizeof(int));
    cdata->read(isCollisionType, numberOfReaction * sizeof(bool));
}

void StoichiometricEquation::Write(DataContainer *cdata)
{
    for ( int irection = 0; irection < numberOfReaction; ++ irection )
    {
        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            cdata->append(&cvp[irection][ispecies], sizeof(int));
        }
    }

    for ( int irection = 0; irection < numberOfReaction; ++ irection )
    {
        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            cdata->append(&cvn[irection][ispecies], sizeof(int));
        }
    }

    for ( int irection = 0; irection < numberOfReaction; ++ irection )
    {
        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            cdata->append(&cpz[irection][ispecies], sizeof(RDouble));
        }
    }

    cdata->append(reactionType, numberOfReaction * sizeof(int));
    cdata->append(isCollisionType, numberOfReaction * sizeof(bool));
}

CurveFitsMethod::CurveFitsMethod()
{
    blotter_ai = 0;
    blotter_bi = 0;
    blotter_ci = 0;
    blotter_di = 0;
    blotter_ei = 0;

    omega_ai = 0;
    omega_bi = 0;
    omega_ci = 0;
    omega_di = 0;

    curvefit_av = 0;
    curvefit_bv = 0; 
    curvefit_cv = 0;

    curvefit_aes = 0;
    curvefit_bes = 0; 
    curvefit_ces = 0;

    numberOfSpecies = 0;
}

CurveFitsMethod::~CurveFitsMethod()
{
    if (numberOfSpecies > 0)
    {
        delete [] blotter_ai;
        delete [] blotter_bi;
        delete [] blotter_ci;
        delete [] blotter_di;
        delete [] blotter_ei;

        delete [] omega_ai;
        delete [] omega_bi;
        delete [] omega_ci;
        delete [] omega_di;

        delete [] curvefit_av;
        delete [] curvefit_bv; 
        delete [] curvefit_cv;

        delete [] curvefit_aes;
        delete [] curvefit_bes;
        delete [] curvefit_ces;

        numberOfSpecies = 0;
    }
}

void CurveFitsMethod::Init(int numberOfSpecies, int nGasModel/* = 0*/)
{
    if (numberOfSpecies <= 0)
    {
        return;
    }
    this->numberOfSpecies = numberOfSpecies;

    blotter_ai = new RDouble[numberOfSpecies];
    blotter_bi = new RDouble[numberOfSpecies];
    blotter_ci = new RDouble[numberOfSpecies];
    blotter_di = new RDouble[numberOfSpecies];
    blotter_ei = new RDouble[numberOfSpecies];

    omega_ai = new RDouble[numberOfSpecies];
    omega_bi = new RDouble[numberOfSpecies];
    omega_ci = new RDouble[numberOfSpecies];
    omega_di = new RDouble[numberOfSpecies];

    curvefit_aes = new RDouble[numberOfSpecies];
    curvefit_bes = new RDouble[numberOfSpecies];
    curvefit_ces = new RDouble[numberOfSpecies];

    //! Set curve fit coefficients of vibrational excitation rate of electron collision for Nitrogen.
    SetNitrogenVibrationalExcitationRate();

    //! Set the curve fit coefficients of effective collision cross-sectional area between the electron and neutral heavy particle.
    SetEffectiveCollisionCrossSectionalAreaCurveFits(nGasModel);
}

//! Set the curve fit coefficients of vibrational excitation rate of electron collision-N2.
void CurveFitsMethod::SetNitrogenVibrationalExcitationRate()
{
    RDouble av[10] = {8.034, 7.924, 7.876, 7.626, 7.326, 4.9, 2.457, 1.119, 0.4681, 0.1837};
    RDouble bv[10] = {-2.227, -2.235, -2.257, -2.334, -2.454, -2.556, -2.702, -2.865, -3.042, -3.223};
    RDouble cv[10] = {2.005, 1.479, 1.054, 0.6499, 0.2049, 7.448e-3, 2.952e-3, 1.133e-3, 4.312e-3, 2.219e-4};

    int nLevel = 10;
    curvefit_av = new RDouble[nLevel];
    curvefit_bv = new RDouble[nLevel];
    curvefit_cv = new RDouble[nLevel];
    for (int i = 0; i < nLevel; i++)
    {
        curvefit_av[i] = av[i];
        curvefit_bv[i] = bv[i];
        curvefit_cv[i] = cv[i];
    }
}

//! Set the curve fit coefficients of effective collision cross-sectional area between the electron and neutral heavy particle.
void CurveFitsMethod::SetEffectiveCollisionCrossSectionalAreaCurveFits(int nGasModel /*= 0*/)
{
    //Curve fitting data for electronic collision area.
    RDouble aes[8] = {1.2e-20, 2.0e-20, 1.0e-19, 5.0e-20, 7.5e-20, 1.0e-20, 1.0e-20, 1.0e-20};
    RDouble bes[8] = {1.7e-24, 6.0e-24, 0.0, 0.0, 5.5e-24, 0.0, 0.0, 0.0};
    RDouble ces[8] = {-2.0e-29, 0.0, 0.0, 0.0, -1.0e-28, 0.0, 0.0, 0.0};

    if (nGasModel == 0) //Earth Gas.
    {
        int nStart = MIN(4, numberOfSpecies);
        for (int i = 0; i < nStart; ++ i)
        {
            curvefit_aes[i] = aes[i];
            curvefit_bes[i] = bes[i];
            curvefit_ces[i] = ces[i];
        }
        for (int i = nStart; i < numberOfSpecies; ++ i)
        {
            curvefit_aes[i] = 0.0;
            curvefit_bes[i] = 0.0;
            curvefit_ces[i] = 0.0;
        }
        if (numberOfSpecies <= 5)
        {
            curvefit_aes[numberOfSpecies - 1] = aes[4];
            curvefit_bes[numberOfSpecies - 1] = bes[4];
            curvefit_ces[numberOfSpecies - 1] = ces[4];
        }
        else
        {
            curvefit_aes[numberOfSpecies - 2] = aes[4];
            curvefit_bes[numberOfSpecies - 2] = bes[4];
            curvefit_ces[numberOfSpecies - 2] = ces[4];
        }
    }
    else
    {
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            curvefit_aes[i] = aes[i];
            curvefit_bes[i] = bes[i];
            curvefit_ces[i] = ces[i];
        }
    }
}

void CurveFitsMethod::SetEffectiveCollisionCrossSectionalAreaCurveFits(vector<string>&namelist)
{
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        if (namelist[s] == "O")
        {
            curvefit_aes[s] = 1.2e-20;
            curvefit_bes[s] = 1.7e-24;
            curvefit_ces[s] = -2.0e-29;
        }
        else if (namelist[s] == "O2")
        {
            curvefit_aes[s] = 2.0e-20;
            curvefit_bes[s] = 6.0e-24;
            curvefit_ces[s] = 0.0;
        }
        else if (namelist[s] == "NO")
        {
            curvefit_aes[s] = 1.0e-19;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
        else if (namelist[s] == "N")
        {
            curvefit_aes[s] = 5.0e-20;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
        else if (namelist[s] == "N2")
        {
            curvefit_aes[s] = 7.5e-20;
            curvefit_bes[s] = 5.5e-24;
            curvefit_ces[s] = -1.0e-28;
        }
        else if (namelist[s] == "C")
        {
            curvefit_aes[s] = 1.0e-20;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
        else if (namelist[s] == "CO")
        {
            curvefit_aes[s] = 1.0e-20;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
        else if (namelist[s] == "CO2")
        {
            curvefit_aes[s] = 1.0e-20;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
        else
        {
            curvefit_aes[s] = 0.0;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
    }
}

void CurveFitsMethod::Read(fstream &file)
{
    string line,word;
    string separator  = " =\t\r\n#$,;\"'";
    string::size_type npos = string::npos;

    //! Obtain the curve fit coefficients of Gupta and Blottner methods for computing the viscosities of species.
    SkipLines(file,2);
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        SkipLines(file,3);
        getline(file,line);

        line = FindNextWord(line,word,separator);
        from_string< RDouble >(blotter_ai[ispecies],word,std::dec);

        line = FindNextWord(line,word,separator);
        from_string< RDouble >(blotter_bi[ispecies],word,std::dec);

        line = FindNextWord(line,word,separator);
        from_string< RDouble >(blotter_ci[ispecies],word,std::dec);

        line = FindNextWord(line,word,separator);
        from_string< RDouble >(blotter_di[ispecies],word,std::dec);

        line = FindNextWord(line,word,separator);
        from_string< RDouble >(blotter_ei[ispecies],word,std::dec);
    }

    //! Obtain the curve fit coefficients of computing the average collision areas for species couples called PI*Omega(2,2).
    SkipLines(file,2);
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        SkipLines(file,3);
        getline(file,line);

        line = FindNextWord(line,word,separator);
        from_string< RDouble >(omega_ai[ispecies],word,std::dec);

        line = FindNextWord(line,word,separator);
        from_string< RDouble >(omega_bi[ispecies],word,std::dec);

        line = FindNextWord(line,word,separator);
        from_string< RDouble >(omega_ci[ispecies],word,std::dec);

        line = FindNextWord(line,word,separator);
        from_string< RDouble >(omega_di[ispecies],word,std::dec);
    }
}

void CurveFitsMethod::GenerateData(vector<string>&namelist)
{
    TransportFitData myCurveFitData;
    TransportFitData *currentData = new TransportFitData();

    myCurveFitData.GenerateFullCurveFittingData();
    myCurveFitData.CreateTransportFitData(namelist, currentData);

    RDouble **viscosityFitData = currentData->GetViscosityCurveFittingData();
    RDouble **collisionAreaFitData = currentData->GetCollisionAreaCurveFittingData();
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        blotter_ai[s] = viscosityFitData[s][0];
        blotter_bi[s] = viscosityFitData[s][1];
        blotter_ci[s] = viscosityFitData[s][2];
        blotter_di[s] = viscosityFitData[s][3];
        blotter_ei[s] = viscosityFitData[s][4];

        omega_ai[s] = collisionAreaFitData[s][0];
        omega_bi[s] = collisionAreaFitData[s][1];
        omega_ci[s] = collisionAreaFitData[s][2];
        omega_di[s] = collisionAreaFitData[s][3];
    }

    SetEffectiveCollisionCrossSectionalAreaCurveFits(namelist);
    delete currentData;
}

void CurveFitsMethod::Read(DataContainer *cdata)
{
    cdata->read(blotter_ai, numberOfSpecies * sizeof(RDouble));
    cdata->read(blotter_bi, numberOfSpecies * sizeof(RDouble));
    cdata->read(blotter_ci, numberOfSpecies * sizeof(RDouble));
    cdata->read(blotter_di, numberOfSpecies * sizeof(RDouble));
    cdata->read(blotter_ei, numberOfSpecies * sizeof(RDouble));

    cdata->read(omega_ai, numberOfSpecies * sizeof(RDouble));
    cdata->read(omega_bi, numberOfSpecies * sizeof(RDouble));
    cdata->read(omega_ci, numberOfSpecies * sizeof(RDouble));
    cdata->read(omega_di, numberOfSpecies * sizeof(RDouble));

    cdata->read(curvefit_av, 10 * sizeof(RDouble));
    cdata->read(curvefit_bv, 10 * sizeof(RDouble));
    cdata->read(curvefit_cv, 10 * sizeof(RDouble));

    cdata->read(curvefit_aes, numberOfSpecies * sizeof(RDouble));
    cdata->read(curvefit_bes, numberOfSpecies * sizeof(RDouble));
    cdata->read(curvefit_ces, numberOfSpecies * sizeof(RDouble));
}

void CurveFitsMethod::Write(DataContainer *cdata)
{
    cdata->append(blotter_ai, numberOfSpecies * sizeof(RDouble));
    cdata->append(blotter_bi, numberOfSpecies * sizeof(RDouble));
    cdata->append(blotter_ci, numberOfSpecies * sizeof(RDouble));
    cdata->append(blotter_di, numberOfSpecies * sizeof(RDouble));
    cdata->append(blotter_ei, numberOfSpecies * sizeof(RDouble));

    cdata->append(omega_ai, numberOfSpecies * sizeof(RDouble));
    cdata->append(omega_bi, numberOfSpecies * sizeof(RDouble));
    cdata->append(omega_ci, numberOfSpecies * sizeof(RDouble));
    cdata->append(omega_di, numberOfSpecies * sizeof(RDouble));

    cdata->append(curvefit_av, 10 * sizeof(RDouble));
    cdata->append(curvefit_bv, 10 * sizeof(RDouble));
    cdata->append(curvefit_cv, 10 * sizeof(RDouble));

    cdata->append(curvefit_aes, numberOfSpecies * sizeof(RDouble));
    cdata->append(curvefit_bes, numberOfSpecies * sizeof(RDouble));
    cdata->append(curvefit_ces, numberOfSpecies * sizeof(RDouble));
}

//void CurveFitsMethod::SetBlotterCurveFittingCoefficients(string *speciesName, int nSpecies)
//{
//    int nIndex = 0;
//    RDouble visCoef[55] = {
//        0.0, 0.0, 0.0203144,  0.4294404, -11.6031403,   //O
//        0.0, 0.0, 0.0449290, -0.0826158, -9.2019475,    //O2
//        0.0, 0.0, 0.0436378, -0.0335511, -9.5767430,    //NO
//        0.0, 0.0, 0.0115572,  0.6031679, -12.432750,    //N
//        0.0, 0.0, 0.0268142,  0.3177838, -11.3155513,   //N2
//        0.0912, -3.3154, 45.1290, -270.4211, 586.2903,  //O+
//        0.0908, -3.3020, 44.9511, -269.3877, 584.4130,  //O2+
//        0.0913, -3.3178,  45.1426, -270.3948, 586.3300, //NO+
//        0.0895, -3.2573,  44.3511, -265.8276, 576.1313, //N+
//        0.0897, -3.2618,  44.4079, -266.1462, 577.1449, //N2+
//        0.0899, -3.2731, 44.5782, -267.2522, 574.4149}; //e-
//
//    for (int s = 0; s < nSpecies; ++ s)
//    {
//        nIndex = GetSpeciesIndex(speciesName[s]);
//        int nStart = 5 * nIndex;
//
//        blotter_ai[s] = visCoef[nStart];
//        blotter_bi[s] = visCoef[nStart + 1];
//        blotter_ci[s] = visCoef[nStart + 2];
//        blotter_di[s] = visCoef[nStart + 3];
//        blotter_ei[s] = visCoef[nStart + 4];
//    }
//}
//
//void CurveFitsMethod::SetGuptaCurveFittingCoefficients(string *speciesName, int nSpecies)
//{
//    int nIndex = 0;
//    RDouble visCoef[55] = {
//        0.0, 0.0, 0.0205, 0.4257, -11.5803, //O
//        0.0, 0.0, 0.0484, -0.1455, -8.9231, //O2
//        0.0, 0.0, 0.0452, -0.0609, -9.4596, //NO
//        0.0, 0.0, 0.0120, 0.5930, -12.3805, //N
//        0.0, 0.0, 0.0203, 0.4329, -11.8153, //N2
//        0.0, 0.0, 0.0000, 2.5000, -32.3606, //O+
//        0.0, 0.0, 0.0000, 2.5000, -32.0148, //O2+
//        0.0, 0.0, 0.0000, 2.5000, -32.0453, //NO+
//        0.0, 0.0, 0.0000, 2.5000, -32.4285, //N+
//        0.0, 0.0, 0.0000, 2.5000, -32.0827, //N2+
//        0.0, 0.0, 0.0000, 2.5000, -37.4475}; //e-
//
//    for (int s = 0; s < nSpecies; ++ s)
//    {
//        nIndex = GetSpeciesIndex(speciesName[s]);
//        int nStart = 5 * nIndex;
//
//        blotter_ai[s] = visCoef[nStart];
//        blotter_bi[s] = visCoef[nStart + 1];
//        blotter_ci[s] = visCoef[nStart + 2];
//        blotter_di[s] = visCoef[nStart + 3];
//        blotter_ei[s] = visCoef[nStart + 4];
//    }
//}

SchmidtNumber::SchmidtNumber()
{
    scl  = 0;
    sct  = 0;
    oscl = 0;
    osct = 0;
    numberOfSpecies = 0;
}

SchmidtNumber::~SchmidtNumber()
{
    delete [] scl ;
    delete [] sct ;
    delete [] oscl;
    delete [] osct;
}

void SchmidtNumber::Init(int numberOfSpecies)
{
    this->numberOfSpecies = numberOfSpecies;

    scl  = new RDouble[numberOfSpecies];
    sct  = new RDouble[numberOfSpecies];
    oscl = new RDouble[numberOfSpecies];
    osct = new RDouble[numberOfSpecies];
}

void SchmidtNumber::ComputeSpeciesSchmidtNumber(int *ionTypeOfSpecies)
{
    //! Calculate the Schmidt number of laminar or turbulent flow of species.

    RDouble sc_l;
    GlobalDataBase::GetData("sc_l", &sc_l, PHDOUBLE, 1);
    RDouble sc_t;
    GlobalDataBase::GetData("sc_t", &sc_t, PHDOUBLE, 1);

    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        if ( ionTypeOfSpecies[ispecies] != 0 && ABS(ionTypeOfSpecies[ispecies]) < 100 )
        {
            scl[ispecies] = half * sc_l;
            sct[ispecies] = half * sc_t;
        }
        else
        {
            scl[ispecies] = 1.0 * sc_l;
            sct[ispecies] = 1.0 * sc_t;
        }

        oscl[ispecies] = 1.0 / scl[ispecies];
        osct[ispecies] = 1.0 / sct[ispecies];
    }

}

void SchmidtNumber::Read(DataContainer *cdata)
{
    cdata->read(scl, numberOfSpecies * sizeof(RDouble));
    cdata->read(sct, numberOfSpecies * sizeof(RDouble));
                                        
    cdata->read(oscl, numberOfSpecies * sizeof(RDouble));
    cdata->read(osct, numberOfSpecies * sizeof(RDouble));
}

void SchmidtNumber::Write(DataContainer *cdata)
{
    cdata->append(scl, numberOfSpecies * sizeof(RDouble));
    cdata->append(sct, numberOfSpecies * sizeof(RDouble));
                                          
    cdata->append(oscl, numberOfSpecies * sizeof(RDouble));
    cdata->append(osct, numberOfSpecies * sizeof(RDouble));
}

MolecularProperty::MolecularProperty()
{
    schmidt_number              = 0;
    molecularWeightDimensional  = 0;
    molecularWeight             = 0;
    oMolecularWeightDimensional = 0;
    oMolecularWeight            = 0;

    characteristicTemperatureOfSpecies = 0;
    mass_fraction_reference = 0;
    name_of_species         = 0;
    ionTypeOfSpecies      = 0;
    collisionCrossSection = 0;

    thermodynamicProperties = 0;
    numberOfSpecies = 0;
}

MolecularProperty::~MolecularProperty()
{
    delete schmidt_number;
    delete [] molecularWeightDimensional ;
    delete [] molecularWeight            ;
    delete [] oMolecularWeightDimensional;
    delete [] oMolecularWeight           ;

    delete [] characteristicTemperatureOfSpecies;
    delete [] mass_fraction_reference;
    delete [] name_of_species        ;
    delete [] ionTypeOfSpecies     ;
    delete [] collisionCrossSection;
    delete [] thermodynamicProperties;
}

void MolecularProperty::Init(int numberOfSpecies)
{
    this->numberOfSpecies = numberOfSpecies;

    molecularWeightDimensional  = new RDouble[numberOfSpecies];
    molecularWeight             = new RDouble[numberOfSpecies];
    oMolecularWeightDimensional = new RDouble[numberOfSpecies];
    oMolecularWeight            = new RDouble[numberOfSpecies];

    characteristicTemperatureOfSpecies = new RDouble[numberOfSpecies];
    name_of_species = new string[numberOfSpecies];

    //! ionTypeOfSpecies is the ion type of species.
    ionTypeOfSpecies        = new int[numberOfSpecies];
    mass_fraction_reference = new RDouble[numberOfSpecies];
    collisionCrossSection = new RDouble[numberOfSpecies];

    schmidt_number = new SchmidtNumber();
    schmidt_number->Init(numberOfSpecies);

    //! To set the thermodynamic properties of species under the multi-temperature model.
    thermodynamicProperties = new Thermo_Param[numberOfSpecies];
}

void MolecularProperty::Read(fstream &file)
{
    string line,word;
    string separator  = " =\t\r\n#$,;\"'";
    string::size_type npos = string::npos;

    SkipLines(file,3);
    //! Obtain the name of species.
    getline(file,line);
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        line = FindNextWord(line,name_of_species[ispecies],separator);
    }

    SkipLines(file,3);
    //! Obtain the ion type of species.
    getline(file,line);
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        line = FindNextWord(line,word,separator);
        from_string< int >(ionTypeOfSpecies[ispecies],word,std::dec);
    }

    SkipLines(file,3);
    //! Obtain the molecular weight of species.
    getline(file,line);
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        line = FindNextWord(line,word,separator);
        from_string< RDouble >(molecularWeightDimensional[ispecies],word,std::dec);
    }

    SkipLines(file,3);
    //! Obtain the mass ratio of species.
    getline(file,line);
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        line = FindNextWord(line,word,separator);
        from_string< RDouble >(mass_fraction_reference[ispecies],word,std::dec);
    }

    SkipLines(file,3);
    //! Obtain the collision cross section of species.
    getline(file,line);
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        line = FindNextWord(line,word,separator);
        from_string< RDouble >(collisionCrossSection[ispecies],word,std::dec);
    }

    SkipLines(file,3);
    //! Obtain the characteristic temperature of species.
    getline(file,line);
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        line = FindNextWord(line,word,separator);
        from_string< RDouble >(characteristicTemperatureOfSpecies[ispecies],word,std::dec);
    }

    schmidt_number->ComputeSpeciesSchmidtNumber(ionTypeOfSpecies);

    //! The thermodynamic properties under the multi-temperature model.
    int nParam = 6;     //! To add 6 parameters including the degeneracy, characteristic temperature etc.
    for (int iparam = 0; iparam < nParam; iparam++) 
    {
        SkipLines(file,3);
        //! Obtain the characteristic temperature of species.
        getline(file,line);
        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            line = FindNextWord(line, word, separator);
            if (iparam == 0)             //! To set the type of species.
            {
                from_string< int >(thermodynamicProperties[ispecies].nType, word, std::dec);
            }else if (iparam == 1)       //! To set the first degeneracy of species.
            {
                from_string< int >(thermodynamicProperties[ispecies].gs0, word, std::dec);
            }else if (iparam == 2)       //! To set the second degeneracy of species.
            {
                from_string< int >(thermodynamicProperties[ispecies].gs1, word, std::dec);
            }else if (iparam == 3)       //! To set the characteristic temperature of species in the vibration mode.
            {
                from_string< RDouble >(thermodynamicProperties[ispecies].Tvs.Tv[0], word, std::dec);
            }else if (iparam == 4)       //! To set the characteristic temperature of species in the electron mode.
            {
                from_string< RDouble >(thermodynamicProperties[ispecies].Tes, word, std::dec);
            }else if (iparam == 5)       //! To set the formation enthalpies of species.
            {
                from_string< RDouble >(thermodynamicProperties[ispecies].hs0, word, std::dec);
            }
        }
    }         //! End of reading the parameters to species.
}
void MolecularProperty::GenerateData(vector<string>&namelist, vector<RDouble> &initMassFractions)
{
    AtmosphereData *baseData = new AtmosphereData();
    AtmosphereData *currentData = new AtmosphereData();

    baseData->GenerateFullPhysicalChemicalData();
    baseData->CreateAtmosphereData(namelist, currentData);

    //Obtain the new data according to the name list.
    int *speciesType = currentData->GetSpeciesType();
    int *speciesCharge = currentData->GetSpeciesCharge();
    int *speciesG0 = currentData->GetSpeciesDegeneracy0();
    int *speciesG1 = currentData->GetSpeciesDegeneracy1();
    RDouble *speciesMass = currentData->GetSpeciesMolecularWeight();
    RDouble *speciesTtr = currentData->GetSpeciesTemperature();
    RDouble *speciesTes = currentData->GetSpeciesElectronTemperature();
    RDouble *speciesSigma = currentData->GetSpeciesCollisionArea();
    RDouble *speciesEnthalpy = currentData->GetSpeciesFormationEnthalpy();
    VibrationModeData *speciesVibrationModeData = currentData->GetSpeciesVibrationalModeData();

    for (int s = 0; s < this->numberOfSpecies; ++ s)
    {
        name_of_species[s] = namelist[s];
        mass_fraction_reference[s] = initMassFractions[s];
        //! Generate the ion type of species.
        ionTypeOfSpecies[s] = speciesCharge[s];
        //! Generate the molecular weight of species.
        molecularWeightDimensional[s] = speciesMass[s];
        //! Generate the collision cross section of species.
        collisionCrossSection[s] = speciesSigma[s];
        //! Generate the characteristic temperature of species.
        characteristicTemperatureOfSpecies[s] = speciesTtr[s];

        thermodynamicProperties[s].nType = speciesType[s];
        thermodynamicProperties[s].gs0 = speciesG0[s];
        thermodynamicProperties[s].gs1 = speciesG1[s];
        thermodynamicProperties[s].Tvs = speciesVibrationModeData[s];
        thermodynamicProperties[s].Tes = speciesTes[s];
        thermodynamicProperties[s].hs0 = speciesEnthalpy[s];
    }
    schmidt_number->ComputeSpeciesSchmidtNumber(ionTypeOfSpecies);

    delete baseData;
    delete currentData;
}

void MolecularProperty::Read(DataContainer *cdata)
{
    schmidt_number->Read(cdata);
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        cdata->read_string(name_of_species[ispecies]);
    }

    cdata->read(ionTypeOfSpecies       , numberOfSpecies * sizeof(int));
    cdata->read(mass_fraction_reference, numberOfSpecies * sizeof(RDouble));

    cdata->read(collisionCrossSection,            numberOfSpecies * sizeof(RDouble));
    cdata->read(characteristicTemperatureOfSpecies,numberOfSpecies * sizeof(RDouble));

    cdata->read(molecularWeightDimensional, numberOfSpecies * sizeof(RDouble));
    cdata->read(molecularWeight           , numberOfSpecies * sizeof(RDouble));

    cdata->read(oMolecularWeightDimensional, numberOfSpecies * sizeof(RDouble));
    cdata->read(oMolecularWeight           , numberOfSpecies * sizeof(RDouble));

    cdata->read(thermodynamicProperties, numberOfSpecies * sizeof(Thermo_Param));
}

void MolecularProperty::Write(DataContainer *cdata)
{
    schmidt_number->Write(cdata);
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        cdata->append_string(name_of_species[ispecies]);
    }

    cdata->append(ionTypeOfSpecies       , numberOfSpecies * sizeof(int));
    cdata->append(mass_fraction_reference, numberOfSpecies * sizeof(RDouble));

    cdata->append(collisionCrossSection             , numberOfSpecies * sizeof(RDouble));
    cdata->append(characteristicTemperatureOfSpecies, numberOfSpecies * sizeof(RDouble));

    cdata->append(molecularWeightDimensional, numberOfSpecies * sizeof(RDouble));
    cdata->append(molecularWeight           , numberOfSpecies * sizeof(RDouble));

    cdata->append(oMolecularWeightDimensional, numberOfSpecies * sizeof(RDouble));
    cdata->append(oMolecularWeight           , numberOfSpecies * sizeof(RDouble));

    cdata->append(thermodynamicProperties, numberOfSpecies * sizeof(Thermo_Param));
}

RDouble * MolecularProperty::GetLaminarSchmidtNumberReciprocal() const
{ 
    return schmidt_number->GetLaminarSchmidtNumberReciprocal();
}

RDouble * MolecularProperty::GetTurbulentSchmidtNumberReciprocal() const
{
    return schmidt_number->GetTurbulentSchmidtNumberReciprocal();
}

MixedGas::MixedGas()
{
    molecularProperty      = 0;
    thermodynamicManager   = 0;
    reactionRate           = 0;
    stoichiometricEquation = 0;
    CurveFits    = 0;
    nGasModel = 0; //nGasModel = 0 indicates the gas of Earth, nGasModel = 1 is for gas of Mars.
    parkVDPower = 0.6;
    nTEnergyModel = 0; //The conventional method to computing temperature energy terms.
    nIsChemicalFreeze = 0;
    nIsSuperCatalytic = 1;
    polynomialCoefT2 = 0;
}

MixedGas::~MixedGas()
{
    if (nchem > 0)
    {
        delete molecularProperty;
        delete thermodynamicManager;
        delete reactionRate;
        delete stoichiometricEquation;
        delete CurveFits;

        DeAllocateWorkingSpace();
        RDouble *prim_inf = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
        delete []prim_inf;
        RDouble *catalyticMassFraction = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("catalyticMassFraction"));
        delete []catalyticMassFraction;
    }
}

RDouble MixedGas::ComputeReferenceTotalEnthalpy()
{
    RDouble totalEnthalpy = 0.0, staticEnthalpy = 0.0;
    RDouble *initMassFraction = this->GetInitMassFraction();
    RDouble initVibrationTemperature = 1.0;
    if (GlobalDataBase::IsExist("freestream_vibration_temperature", PHDOUBLE, 1))
    {
        initVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        initVibrationTemperature = initVibrationTemperature / referenceParameterDimensional.GetTemperature();
    }
    RDouble velocity = referenceParameterDimensional.GetVelocity();
    if (this->ntmodel == 1)
    {
        staticEnthalpy = this->GetMixtureGasEnthalpy(initMassFraction, 1.0);
    }
    else
    {
        staticEnthalpy = this->ComputeMixedGasEnthalpy(1.0, initVibrationTemperature, initVibrationTemperature, initMassFraction);
    }
    totalEnthalpy = (staticEnthalpy + 0.5) * velocity * velocity;

    return totalEnthalpy;
}

void MixedGas::AllocateWorkingSpace()
{
    moleFractions = new RDouble[numberOfSpecies];
    massFractions = new RDouble[numberOfSpecies];

    speciesViscosity = new RDouble[numberOfSpecies];
    speciesWeight = new RDouble[numberOfSpecies];

    speciesCp = new RDouble[numberOfSpecies];
    speciesCv = new RDouble[numberOfSpecies];
    speciesEnthalpy = new RDouble[numberOfSpecies];

    workSpecies = new RDouble[numberOfSpecies];
    speciesEnergyTr = new RDouble[numberOfSpecies];
    speciesEnergyVb = new RDouble[numberOfSpecies];
    speciesEnergyEe = new RDouble[numberOfSpecies];
    speciesConductivity = new RDouble[numberOfSpecies];

    dfsdx   = new RDouble[nl + nchem];
    dfsdy   = new RDouble[nl + nchem];
    dfsdz   = new RDouble[nl + nchem];
    fvis    = new RDouble[nl + nchem + ntmodel - 1];

    maximumSpecies = new RDouble[numberOfSpecies];

    Thermo_Energy_temparay.Cps = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Cvs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Cvtrs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Cvvs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Cves = new RDouble[numberOfSpecies];

    Thermo_Energy_temparay.Hs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Es = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Etrs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Evs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Ees = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Rms = new RDouble[numberOfSpecies];

    Thermo_Energy_temparay.forwardRates = new RDouble [numberOfReaction];
    Thermo_Energy_temparay.backwardRates = new RDouble [numberOfReaction];
    Thermo_Energy_temparay.collisionRates = new RDouble [numberOfReaction];
    Thermo_Energy_temparay.OmegaDerivative = NewPointer2 <RDouble> (numberOfSpecies, ntmodel - 1);

    Thermo_Energy_temparay.partialValues = NewPointer2 <RDouble> (numberOfReaction, nEquation);
    Thermo_Energy_temparay.partialValues2 = NewPointer2 <RDouble> (numberOfReaction, nEquation);

    //Obtain the interpolation coefficients for the polynomial.
    int n = this->thermodynamicManager->GetNumberofPolynomialCoefficient();
    polynomialCoefT2 = NewPointer2<RDouble>(numberOfSpecies, n);
}

void MixedGas::DeAllocateWorkingSpace()
{
    delete []moleFractions;
    delete []massFractions;
    delete []speciesViscosity;
    delete []speciesWeight;

    delete []speciesCp;
    delete []speciesCv;
    delete []speciesEnthalpy;
    delete []workSpecies;

    delete []speciesEnergyTr;
    delete []speciesEnergyVb;
    delete []speciesEnergyEe;
    delete []speciesConductivity;

    delete []dfsdx;
    delete []dfsdy;
    delete []dfsdz;
    delete []fvis;
    
    delete []maximumSpecies;

    delete []Thermo_Energy_temparay.Cps;
    delete []Thermo_Energy_temparay.Cvs;
    delete []Thermo_Energy_temparay.Cvtrs ;
    delete []Thermo_Energy_temparay.Cvvs;
    delete []Thermo_Energy_temparay.Cves;

    delete []Thermo_Energy_temparay.Hs ;
    delete []Thermo_Energy_temparay.Es ;
    delete []Thermo_Energy_temparay.Etrs;
    delete []Thermo_Energy_temparay.Evs;
    delete []Thermo_Energy_temparay.Ees;
    delete []Thermo_Energy_temparay.Rms;

    delete []Thermo_Energy_temparay.forwardRates;
    delete []Thermo_Energy_temparay.backwardRates;
    delete []Thermo_Energy_temparay.collisionRates;
    DelPointer2(Thermo_Energy_temparay.OmegaDerivative);
    DelPointer2(Thermo_Energy_temparay.partialValues);
    DelPointer2(Thermo_Energy_temparay.partialValues2);
    DelPointer2(polynomialCoefT2);
}

void MixedGas::Primitive2Conservative(RDouble *prim, RDouble gama, RDouble Tv, RDouble Te, RDouble *q)
{
    using namespace IDX;
    RDouble em;
    RDouble &density = prim[IR];
    RDouble &um = prim[IU];
    RDouble &vm = prim[IV];
    RDouble &wm = prim[IW];
    RDouble &pressure = prim[IP];
    
    //! Obtain the mass fraction of the last species whose label in the collection is ns-1.
    if ( nchem == 1 )
    {
        NormalizePrimitive(prim);
    }

    //! Obtain the total internal energy Em.
    ComputeInternalEnergy(prim, gama, Tv, Te, em);

    q[IR ] = density;
    q[IRU] = density * um;
    q[IRV] = density * vm;
    q[IRW] = density * wm;
    q[IRE] = density * em;

    for ( int m = nm; m < nl; ++ m )
    {
        q[m] = density * prim[m];
    }
    if (nchem == 1)
    {
        q[nl] = density * prim[nl];
    }
    for (int i = 0; i < ntmodel - 1; i++)       //! multi-temperature model.
    {
        q[nl + nchem + i] = density * prim[nl + nchem + i];
    }
}

void MixedGas::Primitive2ConservativeR(RDouble *prim, RDouble gama, RDouble Tv, RDouble Te, RDouble staticE, RDouble *q)
{
    using namespace IDX;
    RDouble em;
    RDouble &density = prim[IR];
    RDouble &um = prim[IU];
    RDouble &vm = prim[IV];
    RDouble &wm = prim[IW];
    RDouble &pressure = prim[IP];

    RDouble squareVelocity;
    squareVelocity = um * um + vm * vm + wm * wm;
    em = staticE + half * squareVelocity;

    q[IR ] = density;
    q[IRU] = density * um;
    q[IRV] = density * vm;
    q[IRW] = density * wm;
    q[IRE] = density * em;

    for ( int m = nm; m < nEquation; ++ m )
    {
        q[m] = density * prim[m];
    }
}


//void RealGas::Conservative2Primitive(RDouble *q, RDouble gama, RDouble *prim, RDouble *temperature)
//{
//    using namespace IDX;
//    RDouble density, oDensity, um, vm, wm, pressure, totalInternalEnergy, squareVelocity, internalEnergy, deltC,fcs,fmin;
//    int mci;
//    density = q[IR];
//    if (density <= zero)
//    {
//        prim[IR] = density;
//        return;
//    }
//
//    oDensity  = 1.0 / density;
//    um   = q[IU] * oDensity;
//    vm   = q[IV] * oDensity;
//    wm   = q[IW] * oDensity;
//    totalInternalEnergy  = q[IP];    ///den*H
//    squareVelocity   = um * um + vm * vm + wm * wm;
//    internalEnergy  = totalInternalEnergy - half * density * squareVelocity;
//
//    prim[IR] = density;
//    prim[IU] = um;
//    prim[IV] = vm;
//    prim[IW] = wm;
//
//    if (nchem == 1)
//    {  if(nSpeciesLimit==1)
//        {
//            fmin=1.0;
//            for (int m = nm; m < nl; ++ m)
//            {       
//                mci = m-nm;
//                deltC = q[m] * oDensity - prim[m];
//                massFractions[mci] = deltC;
//        
//                if( prim[m] +  deltC < 0  )
//                {
//                    fcs = MIN(MAX(prim[m], 0.1) / (SMALL - deltC), 1.0);
//                    fmin = MIN(fmin, fcs);
//                }
//            }
//             
//            for (int m = nm; m < nl; ++ m)
//            {
//                mci = m-nm;
//                massFractions[mci] = prim[m] + fmin * massFractions[mci];
//                prim[m] = MIN( maximumSpecies[mci], MAX(SMALL,  (1.0 - chemicalRelaxCorf) * prim[m] + chemicalRelaxCorf * massFractions[mci]));
//            }
//            //SetMaximumSpecies
//        }
//        else
//        {
//            for (int m = nm; m < nl; ++ m)
//            {
//                prim[m] =MIN( static_cast<RDouble>(one), MAX(SMALL,  (1.0-chemicalRelaxCorf)* prim[m] + chemicalRelaxCorf * q[m] * oDensity));
//            }
//        }
//        NormalizePrimitive(prim);
//    }
//
//    for (int i = 0; i < ntmodel - 1; ++ i)        //! multi-temperature model.
//    {
//        prim[nl + nchem + i] = MAX(q[nl + nchem + i] * oDensity, 1.0e-30);
//    }
//
//    //GetPressure(prim, gama, internalEnergy, temperature[0], pressure);
//    this->GetTemperatureAndPressure(prim, internalEnergy, temperature, pressure);
//
//    prim[IP] =  MAX(SMALL, pressure);
//}


void MixedGas::Conservative2Primitive(RDouble *q, RDouble gama, RDouble *prim, RDouble *temperature)
{
    using namespace IDX;
    RDouble density, oDensity, um, vm, wm, pressure, totalInternalEnergy, squareVelocity, internalEnergy, deltC,fcs,fmin;
    int mci;
    density = q[IR];
    if (density <= zero)
    {
        prim[IR] = density;
        return;
    }

    oDensity  = 1.0 / density;
    um   = q[IU] * oDensity;
    vm   = q[IV] * oDensity;
    wm   = q[IW] * oDensity;
    totalInternalEnergy  = q[IP];    ///den*H
    squareVelocity   = um * um + vm * vm + wm * wm;
    internalEnergy  = totalInternalEnergy - half * density * squareVelocity;

    prim[IR] = density;
    prim[IU] = um;
    prim[IV] = vm;
    prim[IW] = wm;

    if (nchem == 1)
    {  if(nSpeciesLimit==1)
    {
        fmin=1.0;
        for (int m = nm; m < nl; ++ m)
        {       
            mci = m-nm;
            deltC = q[m] * oDensity - prim[m];
            massFractions[mci] = deltC;

            if( prim[m] +  deltC < 0  )
            {
                fcs = MIN(MAX(prim[m], 0.1) / (SMALL - deltC), 1.0);
                fmin = MIN(fmin, fcs);
            }
        }

        for (int m = nm; m < nl; ++ m)
        {
            mci = m-nm;
            massFractions[mci] = prim[m] + fmin * massFractions[mci];
            prim[m] = MIN( maximumSpecies[mci], MAX(SMALL,  (1.0 - chemicalRelaxCorf) * prim[m] + chemicalRelaxCorf * massFractions[mci]));
        }
        //SetMaximumSpecies
    }
    else
    {
        for (int m = nm; m < nl; ++ m)
        {
            prim[m] =MIN( static_cast<RDouble>(one), MAX(SMALL,  (1.0-chemicalRelaxCorf)* prim[m] + chemicalRelaxCorf * q[m] * oDensity));
        }
    }
    NormalizePrimitive(prim);
    }

    for (int i = 0; i < ntmodel - 1; ++ i)        //! multi-temperature model.
    {
        prim[nl + nchem + i] = MAX(q[nl + nchem + i] * oDensity, 1.0e-30);
    }


    int nIdealState = GlobalDataBase::GetIntParaFromDB("nIdealState");
    if (nIdealState == 1)
    {
        prim[IP] = (gama - one) * internalEnergy;
        //! Compute temperature.
        //RDouble gasConstant = this->GetCoefficientOfStateEquation();
        //temperature[ITT] = prim[IP] / (prim[IR] * gasConstant);
        //temperature[ITV] = temperature[ITT];
        //temperature[ITE] = temperature[ITT];
        GetTemperature(prim, temperature[ITT], temperature[ITV], temperature[ITE]);
    }
    else
    {
        //GetPressure(prim, gama, internalEnergy, temperature[0], pressure);
        this->GetTemperatureAndPressure(prim, internalEnergy, temperature, pressure);
        prim[IP] =  MAX(SMALL, pressure);
    }

}

void MixedGas::Conservative2PrimitiveR(RDouble *q, RDouble gama, RDouble *prim, RDouble *temperature)
{
    using namespace IDX;
    RDouble gasConstant = this->GetUniversalGasConstant();
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();

    RDouble oDensity, squareVelocity, internalEnergy, deltC, fcs, fmin;
    RDouble massReciprocal, ceDivideMe;
    RDouble gasMixtureCv, gasMixtureEnergy, TkLimit;
    int mci, nCount;

    RDouble &density = prim[IR];
    RDouble &um = prim[IU];
    RDouble &vm = prim[IV];
    RDouble &wm = prim[IW];
    RDouble &pressure = prim[IP];

    density = q[IR];    
    oDensity  = 1.0 / density;
    um   = q[IRU] * oDensity;
    vm   = q[IRV] * oDensity;
    wm   = q[IRW] * oDensity;
 
    squareVelocity  = um * um + vm * vm + wm * wm;
    internalEnergy  = q[IRE] * oDensity - half * squareVelocity;

    if(nSpeciesLimit==1)
    {
        fmin=1.0;
        for (int m = nm; m < nl; ++ m)
        {       
            mci = m-nm;
            deltC = q[m] * oDensity - prim[m];
            massFractions[mci] = deltC;
        
            if( prim[m] +  deltC < 0  )
            {
                fcs = MIN(MAX(prim[m], 0.1) / (SMALL - deltC), 1.0);
                fmin = MIN(fmin, fcs);
            }
        }
             
        for (int m = nm; m < nl; ++ m)
        {
            mci = m-nm;
            massFractions[mci] = prim[m] + fmin * massFractions[mci];
            prim[m] = MIN( maximumSpecies[mci], MAX(SMALL,  (1.0 - chemicalRelaxCorf) * prim[m] + chemicalRelaxCorf * massFractions[mci]));
        }
            //SetMaximumSpecies
    }
    else
    {
        for (int m = nm; m < nl; ++ m)
        {
            prim[m] =MIN( static_cast<RDouble>(one), MAX(SMALL,  (1.0-chemicalRelaxCorf)* prim[m] + chemicalRelaxCorf * q[m] * oDensity));
        }
    }
    NormalizePrimitive(prim);

    for (int i = 1; i < ntmodel; ++ i)        //! multi-temperature model.
    {
        prim[nl + i] = MAX(q[nl + i] * oDensity, 1.0e-30);
    }

  /*   internalEnergy = internalEnergy * density;
   this->GetTemperatureAndPressure(prim, internalEnergy, temperature, pressure);
    prim[IP] =  MAX(SMALL, pressure);
    return;*/

    //! Obtain the average molecular weight.
    massReciprocal = 0.0;
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        massFractions[i] = prim[nm + i];
        massReciprocal += massFractions[i] / speciesMass[i];
    }
    
    if (nElectronIndex >= 0)
    {
        ceDivideMe = massFractions[nElectronIndex] / speciesMass[nElectronIndex];
    }
    else
    {
        ceDivideMe = 0.0;
    }

    //compute T and P;

    int nMaxStepTemperature = 5;
    RDouble precisionTemperature = 1.0e-3;

    if(ntmodel == 1)
    {
        RDouble &trTemperature = temperature[mTT];
        nCount = 0;
        while (nCount < nMaxStepTemperature)
        {
            nCount++;
            ComputeMixtureCvAndEnergy(massFractions, trTemperature, gasMixtureCv, gasMixtureEnergy);
           
            TkLimit = (gasMixtureEnergy - internalEnergy) / max(gasMixtureCv,SMALL);

            if (TkLimit > 0.9 * trTemperature)
            {
                TkLimit = 0.9 * trTemperature;
            }
            else if(TkLimit < -9.0 * trTemperature)
            {
                TkLimit = -9.0 * trTemperature;
            }
            trTemperature = MAX(trTemperature - TkLimit, trTemperatureMinNonDimensional);

            if (fabs(TkLimit / trTemperature) <= precisionTemperature)
            {
                break;
            }
        }
        //trTemperature = MAX(trTemperature, SMALL);
        pressure = density * gasConstant * trTemperature * massReciprocal;
    }
    else if(ntmodel == 2)
    {
        RDouble &tveTemperature = temperature[mTV], Eve = prim[nl + 1];
        nCount = 0;
        while (nCount < nMaxStepTemperature)
        {
            nCount++;
            ComputeMixtureCvveAndEve(massFractions, tveTemperature, gasMixtureCv, gasMixtureEnergy);
            TkLimit = (gasMixtureEnergy - Eve) / max(gasMixtureCv,SMALL);
            if (TkLimit > 0.9 * tveTemperature)
            {
                TkLimit = 0.9 * tveTemperature;
            }
            else if(TkLimit < -9.0 * tveTemperature)
            {
                TkLimit = -9.0 * tveTemperature;
            }
            tveTemperature = MAX(tveTemperature - TkLimit, veTemperatureMinNonDimensional);

            if (fabs(TkLimit / tveTemperature) <= precisionTemperature)
            {
                break;
            }
        }

        RDouble Etr = internalEnergy - Eve ,&trTemperature = temperature[mTT];
        ComputeMixtureCvtrAndEtr(massFractions, trTemperature, gasMixtureCv, gasMixtureEnergy);
        /*trTemperature = trTemperature + (Etr - gasMixtureEnergy) / gasMixtureCv;
        if(trTemperature<1e-10)
        {
            trTemperature = MAX(trTemperature, SMALL);
        }*/

        TkLimit = (gasMixtureEnergy - Etr) / gasMixtureCv;
        if (TkLimit > 0.9 * trTemperature)
        {
            TkLimit = 0.9 * trTemperature;
        }
        else if(TkLimit < -9.0 * trTemperature)
        {
            TkLimit = -9.0 * trTemperature;
        }

        trTemperature = MAX(trTemperature - TkLimit, trTemperatureMinNonDimensional);


        trTemperature = MAX(trTemperature, SMALL);
        pressure = density * gasConstant * ( trTemperature * (massReciprocal - ceDivideMe) + tveTemperature * ceDivideMe );
    }
    else
    {
        RDouble &teTemperature = temperature[mTE], Ee = prim[nl + 2];
        nCount = 0;
        while (nCount < nMaxStepTemperature)
        {
            nCount++;
            //Tk1 = teTemperature;
            ComputeMixtureCveAndEe(massFractions, teTemperature, gasMixtureCv, gasMixtureEnergy);

            TkLimit = (gasMixtureEnergy - Ee) / max(gasMixtureCv,SMALL);

            if (TkLimit > 0.9 * teTemperature)
            {
                TkLimit = 0.9 * teTemperature;
            }
            else if(TkLimit < -9.0 * teTemperature)
            {
                TkLimit = -9.0 * teTemperature;
            }
            teTemperature = MAX(teTemperature - TkLimit, veTemperatureMinNonDimensional);

            if (fabs(TkLimit / teTemperature) <= precisionTemperature)
            {
                break;
            }
        }
        //teTemperature = MAX(teTemperature, veTemperatureMinNonDimensional);

        RDouble &tvTemperature = temperature[mTV], Ev = prim[nl + 1];
        nCount = 0;
        while (nCount < nMaxStepTemperature)
        {
            nCount++;

            ComputeMixtureCvvAndEv(massFractions, tvTemperature, gasMixtureCv, gasMixtureEnergy);

            TkLimit = (gasMixtureEnergy - Ev) / max(gasMixtureCv,SMALL);

            if (TkLimit > 0.9 * tvTemperature)
            {
                TkLimit = 0.9 * tvTemperature;
            }
            else if(TkLimit < -9.0 * tvTemperature)
            {
                TkLimit = -9.0 * tvTemperature;
            }
            tvTemperature = MAX(tvTemperature - TkLimit, veTemperatureMinNonDimensional);

            if (fabs(TkLimit / tvTemperature) <= precisionTemperature)
            {
                break;
            }
        }
        //tvTemperature = MAX(tvTemperature, veTemperatureMinNonDimensional);

        RDouble Etr = internalEnergy - Ev - Ee ,&trTemperature = temperature[mTT];
        ComputeMixtureCvtrAndEtr(massFractions, trTemperature, gasMixtureCv, gasMixtureEnergy);
        /*trTemperature = trTemperature + (Etr - gasMixtureEnergy) / gasMixtureCv;
        trTemperature = MAX(trTemperature, SMALL);*/
        /*trTemperature = trTemperature + (Etr - gasMixtureEnergy) / gasMixtureCv;
        if(trTemperature<1e-10)
        {
            trTemperature = MAX(trTemperature, SMALL);
        }*/

        TkLimit = (gasMixtureEnergy - Etr) / gasMixtureCv;
        if (TkLimit > 0.9 * trTemperature)
        {
            TkLimit = 0.9 * trTemperature;
        }
        else if(TkLimit < -9.0 * trTemperature)
        {
            TkLimit = -9.0 * trTemperature;
        }

        trTemperature = MAX(trTemperature - TkLimit, trTemperatureMinNonDimensional);
        pressure = density * gasConstant * ( trTemperature * (massReciprocal - ceDivideMe) + teTemperature * ceDivideMe );
    }
}


RDouble * MixedGas::GetLaminarSchmidtNumberReciprocal() const
{ 
    return molecularProperty->GetLaminarSchmidtNumberReciprocal();
}

RDouble * MixedGas::GetTurbulentSchmidtNumberReciprocal() const
{
    return molecularProperty->GetTurbulentSchmidtNumberReciprocal();
}

//! To obtain the laminar Schmidt number of each component.
void MixedGas::GetLaminarSchmidtNumber(RDouble *Scl) const
{
    if (Scl == NULL)
    {
        return;
    }
    else
    {
        RDouble *pSchmidt = molecularProperty->GetLaminarSchmidtNumberReciprocal();
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            Scl[i] = 1.0 / pSchmidt[i];
        }
    }
}

//! To obtain the turbulent Schmidt number of each component.
void MixedGas::GetTurbulentSchmidtNumber(RDouble *Sct) const
{
    if (Sct == NULL)
    {
        return;
    }
    else
    {
        RDouble *pSchmidt = molecularProperty->GetTurbulentSchmidtNumberReciprocal();
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            Sct[i] = 1.0 / pSchmidt[i];
        }
    }
}

void MixedGas::GetSpecificHeatRatio(RDouble *prim, RDouble &gama)
{
    if (nchem == 0)      //! Perfect gas.
    {
        gama = 1.4;
    }
    else                 //! The chemical non-equilibrium flow.
    {
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            massFractions[i] = prim[nm + i];
        }
        RDouble cp, cv;
        RDouble R = this->GetUniversalGasConstant();
        RDouble M = this->GetMixedGasMolecularWeight(massFractions);
        if (ntmodel >= 2)     //! Two-Temperature and Three-Temperature model.
        {
            cv = this->GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(massFractions);
            cp = cv + R/M;
        }
        else                 //! The default setting is One-Temperature model.
        {
            RDouble Ttr = prim[IDX::IP] * M / (prim[IDX::IR] * R);
            this->ComputeConstantPressureSpecificHeatByMassFraction(massFractions, Ttr, cp);
            cv = cp - R/M;
        }

        gama = cp / cv;
    }
}

void MixedGas::GetSpecificHeatRatioAndTemperatute(RDouble *primitiveVariables, RDouble &gama, RDouble &transRotationTemperature, RDouble &vibrationTemperature, RDouble &electronTemperature)
  {
  
    GetTemperature(primitiveVariables, transRotationTemperature, vibrationTemperature, electronTemperature);

    if (nchem == 0)       //! Perfect gas.
    {
        gama = 1.4;
    }
    else                  //! The chemical non-equilibrium flow.
    {
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            massFractions[i] = primitiveVariables[nm + i];
        }
        RDouble cp, cv;
        RDouble gasConstant = this->GetUniversalGasConstant();
        RDouble mass = this->GetMixedGasMolecularWeight(massFractions);
        //! The NormalizePrimitive() function has been used in the function of GetTemperature().
        if (ntmodel >= 2)        //! Two-Temperature and Three-Temperature model.
        {
            cv = this->GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(massFractions);
            cp = cv + gasConstant / mass;
        }
        else                    //! The default setting is One-Temperature model.
        {
            this->ComputeConstantPressureSpecificHeatByMassFraction(massFractions, transRotationTemperature, cp);
            cv = cp - gasConstant / mass;
        }
        gama = cp / cv;
    }
}

void MixedGas::ReadGasModel()
{
    string gasfile;
    GlobalDataBase::GetData("gasfile", &gasfile, PHSTRING, 1);

    string chemicalModel = "";
    int nSpecies = 0, nReaction = 0;
    if (gasfile == "Gu5")       //! The Gupta 5-species-6-reactions model.
    {
        chemicalModel = "Gupta";
        nSpecies = 5;
        nReaction = 6;
    }
    else if (gasfile == "Gu7")  //! The Gupta 7-species-9-reactions model.
    {
        chemicalModel = "Gupta";
        nSpecies = 7;
        nReaction = 9;
    }
    else if (gasfile == "Gu11") //! The Gupta 11-species-20-reactions model.
    {
        chemicalModel = "Gupta";
        nSpecies = 11;
        nReaction = 20;
    }
    else if (gasfile == "Pa5")  //! The Park 5-species-17-reactions model.
    {
        chemicalModel = "Park";
        nSpecies = 5;
        nReaction = 17;
    }
    else if (gasfile == "Pa7")  //! The Park 7-species-22-reactions model.
    {
        chemicalModel = "Park";
        nSpecies = 7;
        nReaction = 22;
    }
    else if (gasfile == "Pa11") //! The Park 11-species-48-reactions model.
    {
        chemicalModel = "Park";
        nSpecies = 11;
        nReaction = 48;
    }
    else if (gasfile == "DK5")  //! The Dunn-Kang 5-species-11-reactions model.
    {
        chemicalModel = "Dunn-Kang";
        nSpecies = 5;
        nReaction = 11;
    }
    else if (gasfile == "DK7")  //! The Dunn-Kang 7-species-15-reactions model.
    {
        chemicalModel = "Dunn-Kang";
        nSpecies = 7;
        nReaction = 15;
    }
    else if (gasfile == "DK11") //! The Dunn-Kang 11-species-26-reactions model.
    {
        chemicalModel = "Dunn-Kang";
        nSpecies = 11;
        nReaction = 26;
    }
    else if (gasfile == "Mars-Pa8") //! The Park 8-species-12-reactions model for Mars.
    {
        chemicalModel = "Mars-Park";
        nSpecies = 8;
        nReaction = 12;
    }
    else if (gasfile == "Mars-Mc8") //! The McKenzie 8-species-12-reactions model for Mars.
    {
        chemicalModel = "Mars-McKenzie";
        nSpecies = 8;
        nReaction = 12;
    }
    else if (gasfile == "Combustion-12")
    {
        chemicalModel = "CombustionGas"; //The 12-species-20-reactions model of Combustion chamber Gas which reference AEROPH.
        nSpecies = 12;
        nReaction = 20;
    }
    else if (gasfile == "Gas-Mixture")
    {
        chemicalModel = "Gas-Mixture";
        nSpecies = 2;
        nReaction = 1;
    }

    else
    {
        if (gasfile.find("Gu") != gasfile.npos)
        {
            chemicalModel = "Gupta";
        }
        else if (gasfile.find("Pa") != gasfile.npos)
        {
            if (gasfile.find("Mars-Pa") != gasfile.npos)
            {
                chemicalModel = "Mars-Park";
            }
            else
            {
                chemicalModel = "Park";
            }
        }
        else if (gasfile.find("DK") != gasfile.npos)
        {
            chemicalModel = "Dunn-Kang";
        }
        else if (gasfile.find("Mars-Mc") != gasfile.npos)
        {
            chemicalModel = "Mars-McKenzie";
        }
    }

    if (chemicalModel == "")    //! Obtain the gas model from file.
    {
        fstream file;
        file.open(gasfile.c_str(), ios_base::in);

        ReadGasModel(file);

        file.close();
        file.clear();
    }
    else                        //! Generate the gas model in the code.
    {
        GenerateGasModel(chemicalModel);
    }

}

void MixedGas::InitGasModel()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if ( nchem <= 0 ) return;

    if ( myid == server )
    {
        ReadGasModel();
    }

    int number_of_processor = GetNumberOfProcessor();

    if ( number_of_processor <= 1 ) return;

    int tag = 0;
    PH_BcastComposite(ChemCompressData,ChemDecompressData,server,tag);
}

void MixedGas::ComputeReferenceViscosityDimensionalChemical()
{
    RDouble reference_viscosity_dimensional;
    RDouble reference_tempereature = 1.0;

    RDouble *mass_fraction_reference = this->molecularProperty->GetInitialMassFraction();
    ComputeMoleFractionByMassFraction( mass_fraction_reference, moleFractions );

    if (ntmodel > 1)       //! Multi-Temperature Model.
    {
        int nType = 0;     //! 0 for Blottnner curve fit method, 1 for Lennard-Jones method.
        //Obtain the partition function values.
        RDouble *MolecularWeight = this->molecularProperty->GetMolecularWeight();

        if (nType == 0)
        {
            //! Employ Blottnner curve fit method to compute viscosity.
            ComputeSpeciesDimensionalViscosityWithCurveFitMethod(reference_tempereature, reference_tempereature, speciesViscosity);
        }
        else
        {
            //! The reference electron pressure.
            RDouble Pe = 0.0;
            int nIndex = this->nElectronIndex;
            if (nIndex >= 0)
            {
                Pe = 1.0 * mass_fraction_reference[nIndex] * GetUniversalGasConstant() * reference_tempereature / MolecularWeight[nIndex];
            }

            //! Employ Lennard-Jones method to compute viscosity.
            ComputeSpeciesDimensionalViscosityWithLennardJonesMethod(reference_tempereature, Pe, speciesViscosity);
        }

        //! Obtain the partition function values.
        this->GetPartitionFunctionValues(speciesViscosity, MolecularWeight, moleFractions, speciesWeight);

        //! Compute the viscosity of the mixed gas.
        reference_viscosity_dimensional = GetMixedGasViscosityWithWilkeFormula(moleFractions, speciesViscosity, speciesWeight);
    }
    else         //! The default model: one-temperature model.
    { 
        ComputeOneTemperatureModelSpeciesViscosityDimensional( reference_tempereature, speciesViscosity );
        ComputeMixtureCoefficientByWilkeFormula( moleFractions, speciesViscosity, speciesWeight );
        ComputeMixtureByWilkeFormula( moleFractions, speciesViscosity, speciesWeight, reference_viscosity_dimensional );
    }

    referenceParameterDimensional.SetViscosity( reference_viscosity_dimensional );
}

void MixedGas::ComputeReferenceViscosityDimensionalNoChemical()
{
    RDouble t0i  = 288.15;
    RDouble tsi  = 110.4;
    RDouble mu0i = 1.7894E-05;

    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();
    RDouble reference_viscosity_dimensional = (t0i + tsi) / (refDimensionalTemperature + tsi) * pow(refDimensionalTemperature / t0i, 1.5) * mu0i;
    referenceParameterDimensional.SetViscosity(reference_viscosity_dimensional);
}

void MixedGas::ComputeReferenceTsuth()
{
    RDouble tsi  = 110.4;

    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();

    RDouble tsuth = tsi / refDimensionalTemperature;
    GlobalDataBase::UpdateData("tsuth", &tsuth, PHDOUBLE, 1);
}

void MixedGas::ComputeReferenceViscosityDimensional()
{
    if ( nchem > 0 )
    {
        ComputeReferenceViscosityDimensionalChemical();
    }
    else
    {
        ComputeReferenceViscosityDimensionalNoChemical();
    }

    ComputeReferenceTsuth();
}

void MixedGas::ComputeReferenceSpecificHeatRatio()
{
    if ( nchem > 0 )
    {
        ComputeReferenceSpecificHeatRatioWithChemical();
    }
    else
    {
        RDouble refGama = 1.4;
        refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
        referenceParameterDimensional.SetAverageSpecificHeatRatio(refGama);
    }
}

void MixedGas::ComputeReferenceSpecificHeatRatioWithChemical()
{
    RDouble reference_temperature = 1.0;
    RDouble refGama = 1.4;
    RDouble *mass_fraction_reference = this->molecularProperty->GetInitialMassFraction();

    if (ntmodel > 1)       //! multi-temperature model.
    {
        RDouble cvtr = this->GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(mass_fraction_reference);
        //RDouble cvv = this->GetMixedGasDimensionalVibrationSpecificHeatAtConstantVolume(reference_temperature, mass_fraction_reference);
        //RDouble cve = this->GetMixedGasDimensionalElectronSpecificHeatAtConstantVolume(reference_temperature, mass_fraction_reference);

        RDouble cv = cvtr;
        RDouble Mref = referenceParameterDimensional.GetAverageMolecularWeight();
        RDouble cp = cv + rjmk / Mref;
        refGama = cp / cv;
    }
    else                  //! one-temperature model.
    {
        RDouble cp_dimensional, cv_dimensional;
        ComputeSpeciesConstantPressureSpecificHeatDimensional(reference_temperature, speciesCp);
        ComputeMixtureByMassFraction(mass_fraction_reference, speciesCp, cp_dimensional);

        RDouble reference_average_molecular_weight_dimensional   = referenceParameterDimensional.GetAverageMolecularWeight();
        RDouble o_reference_average_molecular_weight_dimensional = 1.0 / reference_average_molecular_weight_dimensional;

        cv_dimensional = cp_dimensional - rjmk * o_reference_average_molecular_weight_dimensional;
        refGama        = cp_dimensional / cv_dimensional;
    }

    referenceParameterDimensional.SetAverageSpecificHeatRatio(refGama);
    GlobalDataBase::UpdateData("refGama", &refGama, PHDOUBLE, 1);
}

void MixedGas::ComputeReferenceGeneralGasConstant()
{
    RDouble reference_average_molecular_weight_dimensional = referenceParameterDimensional.GetAverageMolecularWeight();
    RDouble reference_average_general_gas_constant = rjmk / reference_average_molecular_weight_dimensional;

    referenceParameterDimensional.SetAverageGeneralGasConstant( reference_average_general_gas_constant );
    GlobalDataBase::UpdateData("reference_average_general_gas_constant", &reference_average_general_gas_constant, PHDOUBLE, 1);
}

void MixedGas::SetAirInformationByDataBase()
{
    RDouble height = GlobalDataBase::GetDoubleParaFromDB("height");

    RDouble refDimensionalTemperature;
    RDouble refDimensionalPressure;
    RDouble refDimensionalDensity;
    RDouble refDimensionalSonicSpeed;

    GetAirInfo(height, refDimensionalTemperature, refDimensionalPressure, refDimensionalDensity, refDimensionalSonicSpeed);

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);

    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void MixedGas::SetAirInformationByExpData()
{
    RDouble refDimensionalTemperature;
    RDouble refDimensionalPressure;
    RDouble refDimensionalDensity;
    RDouble refDimensionalSonicSpeed;

    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble T0 = referenceParameterDimensional.GetTemperature();
    RDouble P0 = referenceParameterDimensional.GetPressure();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble R = 287.053;

    RDouble mi = refGama / (refGama - 1.0);
    RDouble coef = 1.0 + (refGama - 1.0) * refMachNumber * refMachNumber / 2.0;

    RDouble T, P, rho, c;
    T = T0 / coef;
    P = P0 / pow( coef, mi );
    rho = P / R / T;
    c = sqrt(refGama * R * T);

    refDimensionalDensity = rho;
    refDimensionalPressure = P;
    refDimensionalTemperature = T;
    refDimensionalSonicSpeed = c;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);

    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void MixedGas::ComputeReferenceGasInformation()
{
    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");

    if (inflowParaType == NONDIMENSIONCONDITION)    //! Obtain the reference density and reference pressure for non-dimension condition.
    {
        ComputeDensityWithReynoldsNumber();
        ComputePressureInGasStateEquation();
    }
    else if ( inflowParaType == FLIGHTCONDITION || inflowParaType == EXPERIMENTCONDITION || inflowParaType == SUBSONICBCCONDITION )      //! Obtain the reference density for non-dimension, and the Re.
    {
        NormalizeAirInformation();
        ComputeReferenceReynoldsNumber();

        if ( inflowParaType == SUBSONICBCCONDITION )
        {
            GetSubsonicBCData();
        }
    }
    else if (inflowParaType == TEMPERATURE_DENSITY)
    {
        RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();
        RDouble refDimensionalDensity = referenceParameterDimensional.GetDensity();
        RDouble gasConstant = referenceParameterDimensional.GetAverageGeneralGasConstant();
        
        RDouble refDimensionalPressure = refDimensionalDensity * gasConstant * refDimensionalTemperature;
        referenceParameterDimensional.SetPressure(refDimensionalPressure);
        GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);

        ComputeReferenceReynoldsNumber();
    }
    else if (inflowParaType == TEMPERATURE_PRESSURE)
    {
        NormalizeAirInformation();
        ComputeReferenceReynoldsNumber();
    }
    else
    {
        TK_Exit::UnexpectedVarValue("inflowParaType", inflowParaType);
    }
}

void MixedGas::ComputeReferenceSoundVelocity()
{
    RDouble gama0 = referenceParameterDimensional.GetAverageSpecificHeatRatio();
    RDouble reference_general_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();

    //! Unit: m/s.
    RDouble refDimensionalSonicSpeed = sqrt(gama0 * reference_general_gas_constant * refDimensionalTemperature);

    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void MixedGas::ComputeReferenceMolecularInformation()
{
    if (nchem > 0)
    {
        ComputeReferenceMolecularInformationWithChemical();
    }
    else
    {
        ComputeReferenceMolecularInformationWithoutChemical();
    }

    RDouble refMass = referenceParameterDimensional.GetAverageMolecularWeight();
    this->electronMassDimension = 5.486e-7;     //! Unit: kg/m3.
    this->electronMass = 5.486e-7 / refMass;    //! non-dimensional value.
}

void MixedGas::ComputeReferenceMolecularInformationWithoutChemical()
{
    //! Unit: kg/mol.
    RDouble reference_average_molecular_weight_dimensional = 28.9644 * 1.0e-3;

    referenceParameterDimensional.SetAverageMolecularWeight( reference_average_molecular_weight_dimensional );
}

void MixedGas::ComputeReferenceMolecularInformationWithChemical()
{
    RDouble *molecularWeight  = molecularProperty->GetMolecularWeight          ();
    RDouble *oMolecularWeight = molecularProperty->GetMolecularWeightReciprocal();

    RDouble *molecularWeightDimensional  = molecularProperty->GetMolecularWeightDimensional();
    RDouble *oMolecularWeightDimensional = molecularProperty->GetMolecularWeightReciprocalDimensional();

    RDouble *mass_fraction_reference = this->molecularProperty->GetInitialMassFraction();

    //! The molecular weight of the species.
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        molecularWeightDimensional[ispecies] *= 1.0e-3;
    }

    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        oMolecularWeightDimensional[ispecies] = 1.0 / molecularWeightDimensional[ispecies];
    }

    RDouble o_reference_average_molecular_weight_dimensional = 0.0;
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        o_reference_average_molecular_weight_dimensional += mass_fraction_reference[ispecies] * oMolecularWeightDimensional[ispecies];
    }

    RDouble reference_average_molecular_weight_dimensional = 1.0 / o_reference_average_molecular_weight_dimensional;

    referenceParameterDimensional.SetAverageMolecularWeight( reference_average_molecular_weight_dimensional );

    //! The non-dimensionalization of the molecular weight of species.
    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        molecularWeight[ispecies] = molecularWeightDimensional[ispecies] * o_reference_average_molecular_weight_dimensional;
    }

    for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
    {
        oMolecularWeight[ispecies] = 1.0 / molecularWeight[ispecies];
    }
}

void MixedGas::ComputeReferenceVelocity()
{
    RDouble refDimensionalSonicSpeed = referenceParameterDimensional.GetSoundVelocity();

    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");

    if (inflowParaType == TEMPERATURE_DENSITY || inflowParaType == TEMPERATURE_PRESSURE)
    {
        RDouble refDimensionalVelocity = referenceParameterDimensional.GetVelocity();
        RDouble refMachNumber = refDimensionalVelocity / refDimensionalSonicSpeed;
        referenceParameterDimensional.SetMachNumber(refMachNumber);
        GlobalDataBase::UpdateData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
    }
    else
    {
        RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();
        RDouble refDimensionalVelocity = refDimensionalSonicSpeed * refMachNumber;
        referenceParameterDimensional.SetVelocity(refDimensionalVelocity);
        GlobalDataBase::UpdateData("refDimensionalVelocity", &refDimensionalVelocity, PHDOUBLE, 1);
    }
}

void MixedGas::ComputeReferenceReynoldsNumber()
{
    RDouble refDimensionalDensity = referenceParameterDimensional.GetDensity();
    RDouble reference_velocity_dimensional  = referenceParameterDimensional.GetVelocity();
    RDouble reference_viscosity_dimensional = referenceParameterDimensional.GetViscosity();

    RDouble refReNumber = refDimensionalDensity * reference_velocity_dimensional / reference_viscosity_dimensional;
    RDouble oreynolds = 1.0 / refReNumber;

    referenceParameterDimensional.SetReynoldsNumber(refReNumber);
    GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("oreynolds", &oreynolds, PHDOUBLE, 1);
}

void MixedGas::ComputeProperReynoldsNumberForGrid()
{
    RDouble refReNumber = referenceParameterDimensional.GetReynoldsNumber();
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    refReNumber *= reynoldsReferenceLengthDimensional;
    RDouble oreynolds = 1.0 / refReNumber;

    referenceParameterDimensional.SetReynoldsNumber(refReNumber);
    GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("oreynolds", &oreynolds, PHDOUBLE, 1);
}

void MixedGas::ComputeOtherProperParameterForGrid()
{
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
    RDouble forceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble forceReferenceLengthSpanWise = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
    RDouble forceReferenceArea = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");

    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    RDouble TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    if (GetDim() == THREE_D)
    {
        forceReferenceArea /= reynoldsReferenceLengthDimensional*reynoldsReferenceLengthDimensional;
    }
    else if (GetDim() == TWO_D)
    {
        forceReferenceArea /= reynoldsReferenceLengthDimensional;
    }
    forceReferenceLength /= reynoldsReferenceLengthDimensional;
    forceReferenceLengthSpanWise /= reynoldsReferenceLengthDimensional;
    TorqueRefX /= reynoldsReferenceLengthDimensional;
    TorqueRefY /= reynoldsReferenceLengthDimensional;
    TorqueRefZ /= reynoldsReferenceLengthDimensional;

    GlobalDataBase::UpdateData("forceReferenceLength", &forceReferenceLength, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("forceReferenceLengthSpanWise", &forceReferenceLengthSpanWise, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("forceReferenceArea", &forceReferenceArea, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);

    int isUnsteady = GlobalDataBase::GetIntParaFromDB( "iunsteady" );
    if ( isUnsteady == 1 )
    {
        PrintToWindow( "Unsteady Informations" , "\n" );

        RDouble physicalTimeStep               = GlobalDataBase::GetDoubleParaFromDB( "physicalTimeStep" );
        RDouble reference_velocity_dimensional = referenceParameterDimensional.GetVelocity();

        RDouble unit = reynoldsReferenceLengthDimensional / reference_velocity_dimensional;

        RDouble physicalTimeStep_dimensional = physicalTimeStep * unit;

        PrintToWindow( "  Physical Time Step       :    ", physicalTimeStep_dimensional , "s", "\n" );
    }
}

void MixedGas::ComputeDensityWithReynoldsNumber()
{
    RDouble refReNumber = referenceParameterDimensional.GetReynoldsNumber();
    RDouble reference_viscosity_dimensional = referenceParameterDimensional.GetViscosity();
    RDouble reference_velocity_dimensional = referenceParameterDimensional.GetVelocity();

    RDouble refDimensionalDensity = refReNumber * reference_viscosity_dimensional / reference_velocity_dimensional;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
}

void MixedGas::ComputePressureInGasStateEquation()
{
    RDouble refDimensionalDensity = referenceParameterDimensional.GetDensity();
    RDouble reference_average_generanl_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();

    RDouble refDimensionalPressure = refDimensionalDensity * reference_average_generanl_gas_constant * refDimensionalTemperature;

    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
}

void MixedGas::InitNumberOfSpecies()
{
    if (nchem == 0)
    {
        numberOfSpecies = 0;
        nl = nm;
        nEquation = nm;
    }
    else
    {
        nl = nm + numberOfSpecies - 1;
        nEquation = nm + numberOfSpecies + ntmodel - 1;
    }
    GlobalDataBase::UpdateData("nl", &nl, PHINT, 1);
    GlobalDataBase::UpdateData("numberOfSpecies", &numberOfSpecies, PHINT, 1);

    if (numberOfSpecies > 0)
    {
        this->nElectronIndex = this->GetSpeciesIndex("e-");
        nNitrogenIndex = this->GetSpeciesIndex("N2");
        if (this->nElectronIndex >= 0)
        {
            this->nSpeciesEquation = numberOfSpecies - 2;
            this->nLastSpeciesIndex = numberOfSpecies - 2;
        }
        else
        {
            this->nSpeciesEquation = numberOfSpecies - 1;
            this->nLastSpeciesIndex = numberOfSpecies - 1;
        }

    }
    else
    {
        this->nSpeciesEquation = 0;
        this->nElectronIndex = -1;
        this->nNitrogenIndex = -1;
        this->nLastSpeciesIndex = -1;
    }

    //! This variable is used to indicate the index of the species with the maximum mass fraction in the array of mixture.
    //for Earth air, Nitrogen is the species with the maximum mass fraction, so it is employed as the variable name.
    GlobalDataBase::UpdateData("nNitrogenIndex", &this->nLastSpeciesIndex, PHINT, 1);
    GlobalDataBase::UpdateData("nElectronIndex", &this->nElectronIndex, PHINT, 1);

    if (nchem > 0)
    {
        AllocateWorkingSpace();
    }
}

RDouble MixedGas::ComputeReferenceReynoldsNumberFromTotalPT( void )
{
    RDouble Re0, refMachNumber;
    RDouble T0 = referenceParameterDimensional.GetTemperature();
    RDouble P0 = referenceParameterDimensional.GetPressure();
    GlobalDataBase::GetData("refMachNumber", &refMachNumber, PHDOUBLE, 1);

    RDouble refGama, c;
    GlobalDataBase::GetData("refGama", &refGama, PHDOUBLE, 1);
    c = referenceParameterDimensional.GetSoundVelocity();
    RDouble R = 287.0;
    RDouble coef, T, P, rho, miuInf;
    RDouble mi = refGama / (refGama - 1.0);
    coef = 1.0 + (refGama - 1.0) * refMachNumber * refMachNumber / 2.0;
    T = T0 / coef;
    P = P0 / pow( coef, mi );
    rho = P / R / T;

    RDouble miu0 = 1.7894e-5;
    RDouble T1 = 288.15, C = 110.4;
    miuInf = pow( T / T1, 1.5 ) * ( T1 + C ) / ( T + C ) * miu0;

    Re0 = rho * refMachNumber * c / miuInf;
    return Re0;
}

void MixedGas::Print()
{
    ostringstream oss;
    string *speciesName = GetNameOfSpecies();
    Thermo_Param *param = molecularProperty->GetThermodynamicProperties();
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        oss<<speciesName[s]<<endl;
        oss<<param[s].nType<<", "<<param[s].Tvs.Tv[0]<<", "<<param[s].Tes<<", "<<param[s].gs0<<", "<<param[s].gs1<<", "<<param[s].hs0<<endl;
    }

    RDouble *blotterA = CurveFits->GetBlotterCurveFittingCoefficientA();
    RDouble *blotterB = CurveFits->GetBlotterCurveFittingCoefficientB();
    RDouble *blotterC = CurveFits->GetBlotterCurveFittingCoefficientC();
    RDouble *blotterD = CurveFits->GetBlotterCurveFittingCoefficientD();
    RDouble *blotterE = CurveFits->GetBlotterCurveFittingCoefficientE();
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        oss<<speciesName[s]<<endl;
        oss<<blotterA[s]<<", "<<blotterB[s]<<", "<<blotterC[s]<<", "<<blotterD[s]<<", "<<blotterE[s]<<endl;
    }

    int nTemp = 5, nCoef = 0;
    nCoef = thermodynamicManager->GetNumberofPolynomialCoefficient();
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        oss<<speciesName[s]<<endl;
        for (int i = 0; i < nTemp; ++ i)
        {
            RDouble *coef = thermodynamicManager->GetPolynomialCoefficient(s, i);
            for (int j = 0; j < nCoef; ++ j)
            {
                oss<<coef[j]<<", ";
            }
            oss<<endl;
        }
    }
    oss<<endl;
    int **forward = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
    int **backward = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
    RDouble **thirdBody = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            oss<<forward[r][s]<<", ";
        }
        oss<<endl;
    }
    oss<<endl;
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            oss<<backward[r][s]<<", ";
        }
        oss<<endl;
    }
    oss<<endl;
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            oss<<thirdBody[r][s]<<", ";
        }
        oss<<endl;
    }
    oss<<endl;

    RDouble *f1 = reactionRate->GetForwardReactionRateCoefficientA();
    RDouble *f2 = reactionRate->GetForwardReactionRateCoefficientB();
    RDouble *f3 = reactionRate->GetForwardReactionRateCoefficientC();
    RDouble *b1 = reactionRate->GetBackwardReactionRateCoefficientA();
    RDouble *b2 = reactionRate->GetBackwardReactionRateCoefficientB();
    RDouble *b3 = reactionRate->GetBackwardReactionRateCoefficientC();
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        oss<<f1[r]<<", "<<f2[r]<<", "<<f3[r]<<", "<<b1[r]<<", "<<b2[r]<<", "<<b3[r]<<endl;
    }
    oss<<endl;

    WriteLogFile(oss);
}

void MixedGas::InitReferenceParameter()
{
    InitGasModel();
    InitNumberOfSpecies();

    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");

    if ( inflowParaType == NONDIMENSIONCONDITION )
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

        referenceParameterDimensional.SetMachNumber(refMachNumber);
        referenceParameterDimensional.SetReynoldsNumber(refReNumber);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    }
    else if ( inflowParaType == FLIGHTCONDITION )
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble height = GlobalDataBase::GetDoubleParaFromDB("height");

        referenceParameterDimensional.SetMachNumber(refMachNumber);

        SetAirInformationByDataBase();
    }
    else if ( inflowParaType == EXPERIMENTCONDITION )
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");

        referenceParameterDimensional.SetMachNumber(refMachNumber);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetPressure(refDimensionalPressure);

        SetAirInformationByExpData();
    }
    else if ( inflowParaType == SUBSONICBCCONDITION )
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");

        referenceParameterDimensional.SetMachNumber(refMachNumber);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetPressure(refDimensionalPressure);
    }
    else if ( inflowParaType == TEMPERATURE_DENSITY ) //where temperature and density are fixed.
    {
        RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");

        referenceParameterDimensional.SetVelocity(refDimensionalVelocity);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetDensity(refDimensionalDensity);
    }
    else if ( inflowParaType == TEMPERATURE_PRESSURE ) //where temperature and pressure are fixed.
    {
        RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");

        referenceParameterDimensional.SetVelocity(refDimensionalVelocity);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetPressure(refDimensionalPressure);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("inflowParaType", inflowParaType);
    }

    RDouble prl = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prt = GlobalDataBase::GetDoubleParaFromDB("prt");

    RDouble oprl = 1.0 / prl;
    GlobalDataBase::UpdateData("oprl", &oprl, PHDOUBLE, 1);
    RDouble oprt = 1.0 / prt;
    GlobalDataBase::UpdateData("oprt", &oprt, PHDOUBLE, 1);

    this->oprl = oprl;
    this->oprt = oprt;
}

void MixedGas::ComputeCoefficientOfStateEquation()
{
    //! No matter whether chemical reactions exist or not.
    RDouble gama0 = referenceParameterDimensional.GetAverageSpecificHeatRatio();
    RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();

    //! The variable is the non-dimensional value of 8.314/M£¬M denotes the non-dimensional molecular weight of the mixture gas.
    coefficientOfStateEquation = 1.0 / (gama0 * refMachNumber * refMachNumber); //! Universal gas constant.
    GlobalDataBase::UpdateData("coefficientOfStateEquation", &coefficientOfStateEquation, PHDOUBLE, 1);
}

void MixedGas::NormalizeAirInformation()
{
    RDouble refDimensionalPressure = referenceParameterDimensional.GetPressure();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();
    RDouble reference_average_generanl_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();

    RDouble refDimensionalDensity = refDimensionalPressure / reference_average_generanl_gas_constant / refDimensionalTemperature;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
}

void MixedGas::InitCommonParameter()
{
    nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    nm    = GlobalDataBase::GetIntParaFromDB("nm");
    ntmodel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    nGasModel = GlobalDataBase::GetIntParaFromDB("nGasModel");
    nTEnergyModel = GlobalDataBase::GetIntParaFromDB("nTEnergyModel");
    parkVDPower = GlobalDataBase::GetDoubleParaFromDB("parkVDPower");
    chemicalRelaxCorf = GlobalDataBase::GetDoubleParaFromDB("chemicalRelaxCorf");
    nSpeciesLimit = GlobalDataBase::GetIntParaFromDB("nSpeciesLimit");
    nEnergyRecycle = GlobalDataBase::GetIntParaFromDB("nEnergyRecycle") * nchem;
    veTemperatureMin = GlobalDataBase::GetDoubleParaFromDB("veTemperatureMin");
    nChemcalSourceModified = GlobalDataBase::GetIntParaFromDB("nChemcalSourceModified");
    nChemcalSourceEsMethod = GlobalDataBase::GetIntParaFromDB("nChemcalSourceEsMethod");

    maxViscous = log(1.0);
    lnMaximum =  log(1.0e36);
    if (GlobalDataBase::IsExist("maxViscous", PHDOUBLE, 1))
    {
        maxViscous = log( GlobalDataBase::GetDoubleParaFromDB("maxViscous")*10.0 );
    }

    trTemperatureMin = 10.0;
    if (GlobalDataBase::IsExist("trTemperatureMin", PHDOUBLE, 1))
    {
        trTemperatureMin = GlobalDataBase::GetDoubleParaFromDB("trTemperatureMin");
    }

    if (GlobalDataBase::IsExist("nIsChemicalFreeze", PHINT, 1))
    {
        nIsChemicalFreeze = GlobalDataBase::GetIntParaFromDB("nIsChemicalFreeze");
    }
    if (GlobalDataBase::IsExist("nIsSuperCatalytic", PHINT, 1))
    {
        nIsSuperCatalytic = GlobalDataBase::GetIntParaFromDB("nIsSuperCatalytic");
        if (nIsSuperCatalytic != 0 && nIsSuperCatalytic != 1)
        {
            //cout<<"The parameter \"nIsSuperCatalytic\" does not equal to zero or one"<<endl;
            nIsSuperCatalytic = 1;
            GlobalDataBase::UpdateData("nIsSuperCatalytic", &nIsSuperCatalytic, PHINT, 1);
        }
    }

    if (parkVDPower < 0.0 || parkVDPower > 1.0)
    {
        parkVDPower = 0.6;
        GlobalDataBase::UpdateData("parkVDPower", &parkVDPower, PHDOUBLE, 1);
    }

    reynoldsReferenceLengthDimensional = 1.0;

    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

    reynoldsReferenceLengthDimensional *= gridScaleFactor;

    GlobalDataBase::UpdateData("reynoldsReferenceLengthDimensional", &reynoldsReferenceLengthDimensional, PHDOUBLE, 1);

    mTT = 0;
    if(ntmodel == 3)
    {
        mTV = 1;
        mTE = 2;
    }
    else if(ntmodel ==2)
    {
        mTV = 1;
        mTE = 1;
    }
    else
    {
        mTV = 0;
        mTE = 0;
    }
}

void MixedGas::InitOtherParameter()
{
    int nolstress, nrokplus;
    nolstress = -1;
    nrokplus  = -1;
    GlobalDataBase::UpdateData("nolstress", &nolstress, PHINT, 1);
    GlobalDataBase::UpdateData("nrokplus" , &nrokplus , PHINT, 1);

    int outnstep = 0;
    GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);
    int innstep = 0;
    GlobalDataBase::UpdateData("innstep", &innstep, PHINT, 1);

    RDouble physicalTime = 0.0;
    GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);

    int nIndexTT = 0;
    int nIndexTV = 0;
    int nIndexTE = 0;
    if (ntmodel == 2)
    {
        nIndexTT = 0;
        nIndexTV = 1;
        nIndexTE = 1;
    }
    else if (ntmodel == 3)
    {
        nIndexTT = 0;
        nIndexTV = 1;
        nIndexTE = 2;
    }
    else
    {
        nIndexTT = 0;
        nIndexTV = 0;
        nIndexTE = 0;
    }
    GlobalDataBase::UpdateData("indexOfTT", &nIndexTT, PHINT, 1);
    GlobalDataBase::UpdateData("indexOfTV", &nIndexTV, PHINT, 1);
    GlobalDataBase::UpdateData("indexOfTE", &nIndexTE, PHINT, 1);
}

void MixedGas::ComputeReferenceParameter()
{
    //! The first two variables is temperature independent, 
    //! depend on molecular weight and mass fraction of species if chemical reactions exist.
    ComputeReferenceMolecularInformation();

    ComputeReferenceGeneralGasConstant();

    //! The following five variables depend on temperature, 
    //! and mass fraction of species if chemical reactions exist.
    ComputeReferenceViscosityDimensional();

    ComputeReferenceSpecificHeatRatio();

    ComputeReferenceSoundVelocity();

    ComputeReferenceVelocity();

    ComputeCoefficientOfStateEquation();
    //! The above seven variables should be computed firstly.

    //! The variables in this function are computed by different methods according to "inflowParaType", 
    //! including density, pressure or reynolds.
    ComputeReferenceGasInformation();

    ComputeReferencePrimitive();
}

void MixedGas::InitGlobalParameterOfNSEquation()
{
    //! Obtain the variable of "reynoldsReferenceLengthDimensional" according to "gridScaleFactor".
    InitCommonParameter();

    InitReferenceParameter();

    ComputeReferenceParameter();

    PrintInflowConditionToWindow();

    ComputeProperReynoldsNumberForGrid();

    //! The five reference variables for the whole aerodynamic coefficients integration, which are from the file of "cfd_para.hypara", should be updated with "reynoldsReferenceLengthDimensional".
    ComputeOtherProperParameterForGrid();

    InitOtherParameter();

    //! Read parameters from file of "boundary_condition.hypara". Notice that the five reference variables for the component aerodynamic coefficients integration should be updated with "reynoldsReferenceLengthDimensional".
    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();

    GlobalBoundaryCondition::SetGlobalBCByGlobalDataBase();

    GlobalBoundaryCondition::SetBCDataBaseByGlobalBC();
}

void MixedGas::ComputeReferencePrimitive()
{
    RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();
    RDouble reference_sound_velocity, reference_pressure;
    
    RDouble massoo_ref = 1.0;

    RDouble reference_density     = 1.0;
    RDouble reference_temperature = 1.0;

    reference_sound_velocity = 1.0 / refMachNumber;
    //! coefficientOfStateEquation = 1.0 / (gama0 * refMachNumber * refMachNumber);
    reference_pressure = coefficientOfStateEquation * 
                         reference_density * 
                         reference_temperature / massoo_ref;

    //! Set the global reference value of energy.
    RDouble refVelocity = referenceParameterDimensional.GetVelocity();
    RDouble refTemperature = referenceParameterDimensional.GetTemperature();
    referenceEnergy = refVelocity * refVelocity;
    referenceSpecificHeat = refTemperature / referenceEnergy;

    veTemperatureMinNonDimensional = veTemperatureMin / refTemperature;
    trTemperatureMinNonDimensional = trTemperatureMin / refTemperature;

    RDouble refTv = reference_temperature;
    RDouble ref_vib_energy = 0.0, ref_ele_energy = 0.0;
    if (nchem == 1)
    {
        //ModifyVibrationEnergy("./test.plt");
        if (ntmodel > 1)
        {
            RDouble *massFraction = molecularProperty->GetInitialMassFraction();
            RDouble ce_div_me = 0.0;
            RDouble Mr = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFraction, ce_div_me);
            RDouble Pe = reference_density * reference_temperature * coefficientOfStateEquation * ce_div_me;

            reference_pressure = reference_density * reference_temperature * coefficientOfStateEquation * Mr + Pe;
            RDouble Tvib = reference_temperature;
            if (GlobalDataBase::IsExist("freestream_vibration_temperature", PHDOUBLE, 1))
            {
                Tvib = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
                Tvib = Tvib / referenceParameterDimensional.GetTemperature();
            }

            //! Compute the vibration energy and electron energy.
            ref_vib_energy = this->GetMixedGasVibrationEnergy(Tvib, massFraction);    //! Tve=Tv=Te=Tinf.
            ref_ele_energy = this->GetMixedGasElectronEnergy(Tvib, massFraction);     //! Tve=Tv=Te=Tinf.
        }
        else       //! One-temperature model.
        {
            reference_pressure = coefficientOfStateEquation * 
                                 reference_density * 
                                 reference_temperature / massoo_ref;
        }
    }

    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    GlobalDataBase::UpdateData("attack", &attack, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("sideslip", &sideslip, PHDOUBLE, 1);

    using namespace IDX;
    RDouble *prim_inf = new RDouble[nl + nchem + ntmodel - 1];
    GlobalDataBase::UpdateDataPtr("prim_inf", prim_inf);

    prim_inf[IR] = reference_density;
    prim_inf[IU] = 1.0 * cos(attack) * cos(sideslip);
    prim_inf[IV] = 1.0 * sin(attack) * cos(sideslip);
    prim_inf[IW] = 1.0 *               sin(sideslip);
    prim_inf[IP] = reference_pressure;

    if (nchem > 0)
    {
        RDouble * mass_fraction_reference = this->molecularProperty->GetInitialMassFraction();

        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            prim_inf[nm + ispecies] = mass_fraction_reference[ispecies];
        }

        //! Set the vibration energy and electron energy.
        if (ntmodel == 2)    //! Tve = Tv = Te = Tinf.
        {
            prim_inf[nl + nchem] = ref_vib_energy + ref_ele_energy;
        }
        else if (ntmodel == 3)
        {
            prim_inf[nl + nchem] = ref_vib_energy;
            prim_inf[nl + nchem + 1] = ref_ele_energy;
        }

        RDouble *catalyticMassFraction = new RDouble[numberOfSpecies];
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            catalyticMassFraction[s] = prim_inf[nm + s];
        }

        if (nIsSuperCatalytic == 1)
        {
            if (this->nGasModel == 0) //Gas of Earth.
            {
                GetEarthFullyCatalyticMassFraction(mass_fraction_reference, catalyticMassFraction);
            }
            else if (this->nGasModel == 1) //Gas of Mars.
            {
                GetMarsFullyCatalyticMassFraction(mass_fraction_reference, catalyticMassFraction);
            }
        }

        SetMaximumSpecies(mass_fraction_reference);

        GlobalDataBase::UpdateDataPtr("catalyticMassFraction", catalyticMassFraction);
    }
}

void MixedGas::CreateAllDataClass(int numberOfSpecies, int numberOfReaction)
{
    molecularProperty = new MolecularProperty();
    molecularProperty->Init(numberOfSpecies);

    reactionRate = new ReactionRate();
    reactionRate->Init(numberOfReaction);

    stoichiometricEquation = new StoichiometricEquation();
    stoichiometricEquation->Init(numberOfReaction, numberOfSpecies);

    CurveFits = new CurveFitsMethod();
    CurveFits->Init(numberOfSpecies);

    thermodynamicManager = new ThermodynamicManager();
    thermodynamicManager->Init(numberOfSpecies);
}

void MixedGas::Read(fstream &file)
{
    molecularProperty->Read(file);
    reactionRate->Read(file);
    thermodynamicManager->Read(file);
    CurveFits->Read(file);
    stoichiometricEquation->Read(file);
}

void MixedGas::ReadGasModel(fstream &file)
{
    string line,word;
    string separator  = " =\t\r\n#$,;\"'";
    string::size_type npos = string::npos;

    SkipLines(file,2);

    getline(file,line);
    line = FindNextWord(line,word,separator);

    from_string< int >(numberOfSpecies,word,std::dec);
    line = FindNextWord(line,word,separator);
    from_string< int >(numberOfReaction,word,std::dec);

    CreateAllDataClass(numberOfSpecies, numberOfReaction);
    Read(file);
}

void MixedGas::GenerateGasModel(string chemicalModel)
{
    int nGasType = 0;
    int Ablation = GlobalDataBase::GetIntParaFromDB("nAblation");
    if (chemicalModel == "Gupta")
    {
        if (Ablation == 0)
        {
            nGasType = 0;
        }
        else
        {
            nGasType = 5;    //! Carbon material ablation.
        }
    }
    else if (chemicalModel == "Park")
    {
        if (Ablation == 0)
        {
            nGasType = 1;
        }
        else
        {
            nGasType = 6;    //! Carbon material ablation.
        }
    }
    else if (chemicalModel == "Dunn-Kang")    //! The Dunn-Kang Model.
    {
        if (Ablation == 0)
        {
            nGasType = 2;
        }
        else
        {
            nGasType = 7;    //! Carbon material ablation.
        }
    }
    else if (chemicalModel == "Mars-Park")    //! The Park Model for Mars.
    {
        nGasType = 3;
    }
    else if (chemicalModel == "Mars-McKenzie")    //! The McKenzie Model for Mars.
    {
        nGasType = 4;
    }
    else if (chemicalModel == "CombustionGas")    //! The Combustion chamber Gas Model.
    {
        nGasType = 8;
    }
    else if (chemicalModel == "Gas-Mixture")    //! The Mixed Gas Model.
    {
        nGasType = 9;
    }

    //Obtain the initial data.
    vector<string> namelist;
    vector<RDouble> initMassFraction;
    this->GetNameList(namelist, initMassFraction);

    //! Generate the chemical reaction model.
    ChemicalReactions myChemModel;
    ChemicalReactions *currentReaction = new ChemicalReactions();
    myChemModel.CreateBasicChemicalReactionModel(nGasType);
    myChemModel.CreateChemicalReactionModel(namelist, currentReaction);

    int ns = 0, nr = 0;
    currentReaction->GetChemicalModelInfo(ns, nr);
    if (ns <= 0 || nr <= 0)
    {
        cout<<"Failed to create the gas model!"<<endl;
        delete currentReaction;
        return;
    }

    //! The number of species and reaction are imported from the database.
    numberOfSpecies = ns;
    numberOfReaction = nr;
    //! Allocate memories.
    CreateAllDataClass(numberOfSpecies, numberOfReaction);

    //-------------------------------------------------------------------
    //! Generate the physical and chemical data of species.
    molecularProperty->GenerateData(namelist, initMassFraction);

    //! Generate the curve fitting data of chemical reaction rate.
    RDouble **reactionCoef = currentReaction->GetReactionRateCurveFittingData();
    reactionRate->Create(reactionCoef);

    //! Generate the curve fitting data of chemkin model for computing the specific heat and enthalpy in the one-temperature model.
    thermodynamicManager->GenerateData(namelist);

    //! Generate the curve fitting data of transport coefficients, such as viscosity.
    CurveFits->GenerateData(namelist);

    //! Generate the chemical reaction equations.
    int **reactionFlag = currentReaction->GetReactionType();
    int **forwardCoef = currentReaction->GetForwardStoichiometricCoefficient();
    int **backwardCoef = currentReaction->GetBackwardStoichiometricCoefficient();
    RDouble **collisionCoef = currentReaction->GetThirdBodyCollisionCoefficient();
    stoichiometricEquation->Create(forwardCoef, backwardCoef, collisionCoef, reactionFlag);

    delete currentReaction;
}

void MixedGas::GetNameList(vector<string> &namelist, vector<RDouble> &initMassFraction)
{
    string speciesName[30], line1, line2;
    RDouble massFraction[30];
    GlobalDataBase::GetData("speciesName", &line1, PHSTRING, 1);
    GlobalDataBase::GetData("initMassFraction", &line2, PHSTRING, 1);
    string word;
    string separator  = " =\t\r\n#$,;\"'";
    int nSpecies = 0;
    for (int s = 0; s < 30; ++ s)
    {
        line1 = FindNextWord(line1, word, separator);
        if (word == "")
        {
            break;
        }
        ++ nSpecies;
        speciesName[s] = word;

        line2 = FindNextWord(line2, word, separator);
        from_string<RDouble>(massFraction[s], word, std::dec);
    }
    namelist.clear();
    initMassFraction.clear();
    for (int s = 0; s < nSpecies; ++ s)
    {
        namelist.push_back(speciesName[s]);
        initMassFraction.push_back(massFraction[s]);
    }
}

void MixedGas::Read(DataContainer *cdata)
{
    molecularProperty->Read(cdata);
    reactionRate->Read(cdata);
    stoichiometricEquation->Read(cdata);
    CurveFits->Read(cdata);
    thermodynamicManager->Read(cdata);
}

void MixedGas::Write(DataContainer *cdata)
{
    molecularProperty->Write(cdata);
    reactionRate->Write(cdata);
    stoichiometricEquation->Write(cdata);
    CurveFits->Write(cdata);
    thermodynamicManager->Write(cdata);
}

void MixedGas::CompressData(DataContainer *&cdata)
{
    cdata->append(&numberOfSpecies, sizeof(int));
    cdata->append(&numberOfReaction, sizeof(int));

    Write(cdata);
}

void MixedGas::DecompressData(DataContainer *cdata)
{
    cdata->moveto_begin();

    cdata->read(&numberOfSpecies, sizeof(int));
    cdata->read(&numberOfReaction, sizeof(int));

    CreateAllDataClass(numberOfSpecies, numberOfReaction);
    Read(cdata);
}

void MixedGas::GetTemperatureRangeIndex(const RDouble &temperature_dimensional, int &indexOfTemperatureRange)
{
    this->thermodynamicManager->GetTemperatureRangeIndex(temperature_dimensional, indexOfTemperatureRange);
}

void MixedGas::ComputeSpeciesConstantPressureSpecificHeatDimensional(RDouble temperature, RDouble *cpSpeciesDimensional)
{
    int ntype = 1;    //! 0 stands for applying the old method, 1 for the new method.

    if (ntype == 0)   //! The primary algorithm.
    {
        RDouble t1, t2, t3, t4;

        RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    
        t1 = temperature * referenceTemperature;
        t2 = t1 * t1;
        t3 = t1 * t2;
        t4 = t1 * t3;

        int itemperature;
        GetTemperatureRangeIndex(t1, itemperature);

        RDouble *oMolecularWeightDimensional = molecularProperty->GetMolecularWeightReciprocalDimensional();

        for ( int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies )
        {
            RDouble *coef = this->thermodynamicManager->GetPolynomialCoefficient(ispecies, itemperature);

            cpSpeciesDimensional[ispecies] = coef[0]
                                           + coef[1] * t1
                                           + coef[2] * t2
                                           + coef[3] * t3
                                           + coef[4] * t4;
            cpSpeciesDimensional[ispecies] *= rjmk * oMolecularWeightDimensional[ispecies];
        }
        return;
    }
    else      //! The new method that modified the coefficients in the temperature boundary range.
    {
        RDouble gasConstant = rjmk;      //! The universal gas constant equals to 8.31434 J/mol-K.
        //! The following variables are dimensional values.
         RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
         RDouble *speciesMass = this->molecularProperty->GetMolecularWeightDimensional(); //The mole weight kg/mol.

        int n = this->thermodynamicManager->GetNumberofPolynomialCoefficient();
        if (n < 1)
        {
            return;
        }

        //! Allocate the memory space.
        RDouble coef[7];
        RDouble T1, T2, T3, T4;
        T1 = temperature * referenceTemperature;
        T1 = MAX(T1, 50.0);
        T2 = T1 * T1;
        T3 = T1 * T2;
        T4 = T1 * T3;

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            //! Apply the polynomial by linearly averaging the coefficients on the temperature boundary.
            GetLinearlyAveragedPolynomialCoefficients(ispecies,T1,coef);

            cpSpeciesDimensional[ispecies] = coef[0] + coef[1] * T1 + coef[2] * T2 + coef[3] * T3 + coef[4] * T4;
            cpSpeciesDimensional[ispecies] *= gasConstant / speciesMass[ispecies];
        }
    }
}

RDouble MixedGas::GetMixedGasMolecularWeight(RDouble *massFraction)
{
    RDouble mixedMass = 0.0;
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        mixedMass += massFraction[i] / speciesMass[i];
    }
    mixedMass = 1.0 / mixedMass;

    //! Return the non-dimensional molecular weight of the mixture gas.
    return mixedMass;
}

RDouble MixedGas::GetMixedGasMolecularWeightReciprocalWithoutElectron(RDouble *massFraction, RDouble &ceDivideMe)
{
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();
    //! Obtain the average molecular weight.
    RDouble mixedMass = GetMixedGasMolecularWeight(massFraction);
    ceDivideMe = 0.0;
    if (nElectronIndex >= 0)
    {
        ceDivideMe = massFraction[nElectronIndex] / speciesMass[nElectronIndex];
    }

    return (1.0 / mixedMass - ceDivideMe);
}

RDouble MixedGas::GetUniversalGasConstant()
{
    //! To obtain the non-dimensional universal gas constant called R = 8.314J/(mol*K).
    //RDouble gasConstant = rjmk;//!unit£ºJ/(mol*K)
    //RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    //RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    //RDouble referenceMass = this->referenceParameterDimensional.GetAverageMolecularWeight();
    //gasConstant = gasConstant * referenceTemperature / (referenceMass * referenceVelocity * referenceVelocity);

    //// Return the non-dimensional universal gas constant.
    //return gasConstant;
    return coefficientOfStateEquation;
}

RDouble MixedGas::GetMixtureGasEnthalpy(RDouble *massFraction, RDouble temperature)
{
    //! Obtain the enthalpy of each species.
    this->ComputeSpeciesEnthalpy(temperature, speciesEnthalpy);
    //! Integral.
    RDouble mixedGasEnthalpy = 0.0;
    this->ComputeMixtureByMassFraction(massFraction, speciesEnthalpy, mixedGasEnthalpy);

    return mixedGasEnthalpy;
}

void MixedGas::ComputeConstantVolumeSpecificHeatByMassFraction(RDouble *massFraction, RDouble temperature, RDouble &mixedGasCv)
{
    RDouble mixedGasCp;       //! Specific heat at constant pressure.
    this->ComputeConstantPressureSpecificHeatByMassFraction(massFraction, temperature, mixedGasCp);

    //! Obtain universal gas constant.
    RDouble gasConstant = this->GetUniversalGasConstant();
    //! Obtain molecular weight of the mixture gas.
    RDouble mixedGasMass = this->GetMixedGasMolecularWeight(massFraction);

    //! Compute the constant volume specific heat of the mixture gas.
    mixedGasCv = mixedGasCp - gasConstant / mixedGasMass;
}

void MixedGas::ComputeSpeciesConstantVolumeSpecificHeat(RDouble temperature, RDouble *cvSpecies)
{
    //! To obtain the non-dimensional constant volume specific heat called cv of each species by fitting formula.
    //! temperature is the non-dimensional temperature of the mixture gas.
    //! cvSpecies is an array saving the non-dimensional constant volume specific heat of each species.

    //! Obtain the dimensional value of specific heat at constant volume.
    ComputeSpeciesConstantVolumeSpecificHeatDimensional(temperature,cvSpecies);

    //! Obtain the reference values.
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceValue = referenceTemperature / (referenceVelocity * referenceVelocity);

    //! Translate the dimensional value to the non-dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        cvSpecies[i] = cvSpecies[i] * referenceValue;
    }
}

void MixedGas::ComputeSpeciesConstantVolumeSpecificHeatDimensional(RDouble temperature, RDouble *cvSpeciesDimensional)
{
    //! To obtain the dimensional constant volume specific heat called cv[J/(kg*K)] of each species by fitting formula.
    //! temperature is the non-dimensional temperature of the mixture gas.
    //! cvSpeciesDimensional is an array saving the dimensional constant volume specific heat of each species.
    //! The following variables are the non-dimensional values.

    //! Obtain the dimensional value of constant pressure specific heat of each species.
    this->ComputeSpeciesConstantPressureSpecificHeatDimensional(temperature,cvSpeciesDimensional);

    RDouble gasConstant = rjmk;// dimensional universal gas constant.
    //! Obtain molecular weight of each species.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();

    //! Compute the specific heat at constant volume.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        cvSpeciesDimensional[i] = cvSpeciesDimensional[i] - gasConstant / speciesMass[i];
    }
}

void MixedGas::ComputeConstantPressureSpecificHeatByMassFraction(RDouble *massFraction, RDouble temperature, RDouble &mixedcp)
{
    ComputeSpeciesConstantPressureSpecificHeat(temperature, speciesCp);
    ComputeMixtureByMassFraction(massFraction, speciesCp, mixedcp);
}

void MixedGas::ComputeSpeciesConstantPressureSpecificHeat(RDouble temperature, RDouble *cpSpecies)
{
    //! Obtain the dimensional value of species specific heat at constant pressure.
    this->ComputeSpeciesConstantPressureSpecificHeatDimensional(temperature, cpSpecies);

    //! Change to the non-dimensional values.
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        cpSpecies[ispecies] = cpSpecies[ispecies] * referenceSpecificHeat;
    }
}

void MixedGas::ComputeSpeciesEnthalpy(RDouble temperature, RDouble *speciesEnthaly)
{
    //! Obtain the dimensional values of species enthalpies.
    this->ComputeSpeciesEnthalpyDimensional(temperature, speciesEnthaly);

    //! Change to the non-dimensional values.
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        speciesEnthaly[s] = speciesEnthaly[s] / referenceEnergy;
    }
}

void MixedGas::ComputeSpeciesEnthalpyDimensional(RDouble temperature, RDouble *speciesEnthaly)
{
    int ntype = 1;     //! 0 stands for applying the old method, 1 for the new method.

    if (ntype == 0)    //! The primary algorithm.
    {
        RDouble t1, t2, t3, t4, ot1;

        RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();

        t1 = temperature * referenceTemperature;
        t2 = t1 * t1;
        t3 = t1 * t2;
        t4 = t1 * t3;
        ot1 = 1.0 / t1;

        int itemperature;
        GetTemperatureRangeIndex(t1, itemperature);

        RDouble *oMolecularWeightDimensional = molecularProperty->GetMolecularWeightReciprocalDimensional();

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            RDouble *coef = this->thermodynamicManager->GetPolynomialCoefficient(ispecies, itemperature);
            speciesEnthaly[ispecies] = coef[0]            +
                                       coef[1] * t1 / 2.0 +
                                       coef[2] * t2 / 3.0 +
                                       coef[3] * t3 / 4.0 +
                                       coef[4] * t4 / 5.0 +
                                       coef[5] * ot1;
            speciesEnthaly[ispecies] *= rjmk * oMolecularWeightDimensional[ispecies] * t1;
        }
    }
    else      //! The new method that modified the coefficients in the temperature boundary range.
    {
        RDouble gasConstant = rjmk;     //! The universal gas constant equals to 8.31434 J/mol-K.
        //! The following variables are dimensional values.
        RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
        RDouble *speciesMass = this->molecularProperty->GetMolecularWeightDimensional();      //! The mole weight kg/mol.

        int n = this->thermodynamicManager->GetNumberofPolynomialCoefficient();
        if (n < 1)
        {
            return;
        }

        RDouble coef[7];
        RDouble T1, T2, T3, T4;
        T1 = temperature * referenceTemperature;
        T1 = MAX(T1, 50.0);
        T2 = T1 * T1;
        T3 = T1 * T2;
        T4 = T1 * T3;

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            //! Apply the polynomial by linearly averaging the coefficients on the temperature boundary.
            GetLinearlyAveragedPolynomialCoefficients(ispecies, T1, coef);

            speciesEnthaly[ispecies] = coef[0] + coef[1] * T1/2.0 + coef[2] * T2/3.0 + coef[3] * T3/4.0 + coef[4] * T4/5.0 + coef[5]/T1;
            speciesEnthaly[ispecies] *= T1 * gasConstant / speciesMass[ispecies];
        }

    }
}

void MixedGas::ComputeSpeciesEnergy(RDouble temperature, RDouble *speciesEnergy)
{
    //! Obtain the dimensional values of species energies.
    this->ComputeSpeciesEnergyDimensional(temperature, speciesEnergy);

    //! Change to the non-dimensional values.
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        speciesEnergy[s] = speciesEnergy[s] / referenceEnergy;
    }
}

void MixedGas::ComputeSpeciesEnergyDimensional(RDouble temperature, RDouble *speciesEnergy)
{
    RDouble gasConstant = rjmk;
    RDouble referenceTemperature = referenceParameterDimensional.GetTemperature();
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();

    RDouble speciesEnthalpy[MAX_SPECIES_NUM];
    ComputeSpeciesEnthalpyDimensional(temperature, speciesEnthalpy);
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        speciesEnergy[s] = speciesEnthalpy[s] - gasConstant * temperature * referenceTemperature / speciesMass[s];
    }
}

void MixedGas::ComputeMixtureByMassFraction(RDouble *massFraction, RDouble *mixtureOfSpecies, RDouble &mixture)
{
    mixture = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        mixture += massFraction[ispecies] * mixtureOfSpecies[ispecies];
    }
}

void MixedGas::ComputeMixtureByPrimitive(RDouble *prim, RDouble *mixtureOfSpecies, RDouble &mixture)
{
    mixture = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        mixture += prim[nm + ispecies] * mixtureOfSpecies[ispecies];
    }
}

void MixedGas::ComputeMoleFractionByMassFraction(RDouble *massFractionOfSpecies, RDouble *moleFractionOfSpecies)
{
    RDouble *oMolecularWeight = molecularProperty->GetMolecularWeightReciprocal();     //! The reciprocal value of mass fraction.

    RDouble xiSum, oxiSum;

    xiSum = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        xiSum += massFractionOfSpecies[ispecies] * oMolecularWeight[ispecies];
    }
    oxiSum = 1.0 / xiSum;

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
       moleFractionOfSpecies[ispecies] = massFractionOfSpecies[ispecies] * oMolecularWeight[ispecies] * oxiSum;
    }
}

void MixedGas::ComputeMixtureCoefficientByWilkeFormula(RDouble *moleFractionOfSpecies, RDouble *chem_var, RDouble *chem_phi)
{
//! Arguments: moleFractionOfSpecies saves the mole fraction of the species,chem_var stores the viscosity of the given species,
//! return the value phi stored in the variable chem_phi.

    RDouble *molecularWeight = molecularProperty->GetMolecularWeight();

    RDouble tmp, tmp1, tmp2;

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        chem_phi[ispecies] = 0.0;
        for (int jspecies = 0; jspecies < numberOfSpecies; ++ jspecies)
        {
            tmp = sqrt(molecularWeight[jspecies] / molecularWeight[ispecies]);
            tmp1 = 1.0 + sqrt( chem_var[ispecies] / chem_var[jspecies] ) * sqrt(tmp);
            tmp2 = sqrt( 8.0 * ( 1.0 + molecularWeight[ispecies] / molecularWeight[jspecies] ) );
            chem_phi[ispecies] += moleFractionOfSpecies[jspecies] * SQR(tmp1) / tmp2;
        }
    }
}

void MixedGas::ComputeMixtureByWilkeFormula(RDouble *moleFractionOfSpecies, RDouble *mixtureOfSpecies, RDouble *phiOfSpecies, RDouble &mixture)
{
    mixture = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        mixture += moleFractionOfSpecies[ispecies] * mixtureOfSpecies[ispecies] / phiOfSpecies[ispecies];
    }
}

void MixedGas::ComputeTemperature(RDouble *primitiveVariables, RDouble &temperature)
{
    //! To compute the non-dimensional temperature of the mixture gas in accordance with primitive variables.
    using namespace IDX;
    RDouble massReciprocal, density, pressure;
    NormalizePrimitive(primitiveVariables);
    ComputeMolecularWeightReciprocal(primitiveVariables, massReciprocal);
    density  = primitiveVariables[IR];
    pressure = primitiveVariables[IP];
    
    temperature = pressure / (density * coefficientOfStateEquation * massReciprocal);
}

void MixedGas::GetTemperature(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &vibrationTemperature, RDouble &electronTemperature)
{
    using namespace IDX;
    RDouble density, pressure, gasConstant, Tk = 0.0;
    density = primitiveVariables[IR];
    pressure = primitiveVariables[IP];
    gasConstant = GetUniversalGasConstant();
    //int nChemicalFlowStep = GlobalDataBase::GetIntParaFromDB("nChemicalFlowStep");
    //int nOuterStep = GlobalDataBase::GetIntParaFromDB("outnstep");

    if (nchem == 0)     //! Perfect gas.
    {
        transRotationTemperature = pressure / (density * gasConstant);
    }
    else                //! Chemical non-equilibrium flow.
    {
        //! Compute the mass fraction of the last species.
        NormalizePrimitive(primitiveVariables);

        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            massFractions[i] = primitiveVariables[nm + i];
        }

        RDouble massReciprocal = 0.0, ceDivideMe = 0.0;
        massReciprocal = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFractions, ceDivideMe);

        RDouble  electronPressure;
        if( ntmodel == 1 ) //! The default setting is One-Temperature model.
        {
            transRotationTemperature = pressure / (density * gasConstant * (massReciprocal + ceDivideMe));
        }    
        else if (ntmodel == 2)      //! Two-Temperature model.
        {
            RDouble Eve = primitiveVariables[nl + nchem];
            Tk = vibrationTemperature;
            RDouble Tve = ComputeVibrationAndElectronTemperature(massFractions, Eve, Tk);

            //! Compute pressure of electron.
            electronPressure = density * gasConstant * Tve * ceDivideMe;
            
            transRotationTemperature = (pressure - electronPressure) / (density * gasConstant * massReciprocal);
            vibrationTemperature = Tve;
        }
        else  //! Three-Temperature model.
        {
            RDouble Ev = primitiveVariables[nl + nchem];
            RDouble Ee = primitiveVariables[nl + nchem + 1];
            Tk = vibrationTemperature;
            RDouble Tv = this->ComputeVibrationTemperature(massFractions, Ev, Tk);
            Tk = electronTemperature;
            RDouble Te = this->ComputeElectronTemperature(massFractions, Ee, Tk);
            //! Compute pressure of electron.
            electronPressure = density * gasConstant * Te * ceDivideMe;
            
            transRotationTemperature = (pressure - electronPressure) / (density * gasConstant * massReciprocal);
            vibrationTemperature = Tv;
            electronTemperature = Te;
        }
    }
}

void MixedGas::GetTemperatureR(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &vibrationTemperature, RDouble &electronTemperature)
{
    using namespace IDX;
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble density, pressure, gasConstant;
    density = primitiveVariables[IR];
    pressure = primitiveVariables[IP];
    gasConstant = GetUniversalGasConstant();
    //int nChemicalFlowStep = GlobalDataBase::GetIntParaFromDB("nChemicalFlowStep");
    //int nOuterStep = GlobalDataBase::GetIntParaFromDB("outnstep");

   // Compute the mass fraction of the last species.
    NormalizePrimitive(primitiveVariables);

    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVariables[nm + i];
    }

    RDouble massReciprocal = 0.0, ceDivideMe = 0.0;
    massReciprocal = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFractions, ceDivideMe);

    int nMaxStepTemperature = 5;
    int nCount;
    RDouble gasMixtureCv, TkLimit, gasMixtureEnergy, precisionTemperature = 1.0e-2;

    if(ntmodel == 1)
    {
        transRotationTemperature = pressure / (density * gasConstant * (massReciprocal + ceDivideMe));
    }
    else if(ntmodel == 2)
    {
        RDouble &tveTemperature = vibrationTemperature, Eve = primitiveVariables[nl + nchem];
        nCount = 0;
        while (nCount < nMaxStepTemperature)
        {
            nCount++;
            //Tk1 = tveTemperature;
            ComputeMixtureCvveAndEve(massFractions, tveTemperature, gasMixtureCv, gasMixtureEnergy);
            TkLimit = (gasMixtureEnergy - Eve) / max(gasMixtureCv,SMALL);

            if (TkLimit > 0.9 * tveTemperature)
            {
                TkLimit = 0.9 * tveTemperature;
            }
            else if(TkLimit < -9.0 * tveTemperature)
            {
                TkLimit = -9.0 * tveTemperature;
            }
            tveTemperature = MAX(tveTemperature - TkLimit, veTemperatureMinNonDimensional);

            if (fabs(TkLimit / tveTemperature) <= precisionTemperature)
            {
                break;
            }
        }

        transRotationTemperature = (pressure - density * gasConstant * vibrationTemperature * ceDivideMe) / (density * gasConstant * massReciprocal);
        transRotationTemperature = MAX(transRotationTemperature, SMALL);

    }
    else
    {
        RDouble &teTemperature = electronTemperature, Ee = primitiveVariables[nl + 2];
        nCount = 0;
        while (nCount < nMaxStepTemperature)
        {
            nCount++;
            ComputeMixtureCveAndEe(massFractions, teTemperature, gasMixtureCv, gasMixtureEnergy);

            TkLimit = (gasMixtureEnergy - Ee) / max(gasMixtureCv,SMALL);
            if (TkLimit > 0.9 * teTemperature)
            {
                TkLimit = 0.9 * teTemperature;
            }
            else if(TkLimit < -9.0 * teTemperature)
            {
                TkLimit = -9.0 * teTemperature;
            }
            teTemperature = MAX(teTemperature - TkLimit, veTemperatureMinNonDimensional);

            if (fabs(TkLimit / teTemperature) <= precisionTemperature)
            {
                break;
            }
        }
        //teTemperature = MAX(teTemperature, veTemperatureMinNonDimensional);

        RDouble &tvTemperature = vibrationTemperature, Ev = primitiveVariables[nl + 1];
        nCount = 0;
        while (nCount < nMaxStepTemperature)
        {
            nCount++;
            ComputeMixtureCvvAndEv(massFractions, tvTemperature, gasMixtureCv, gasMixtureEnergy);

            TkLimit = (gasMixtureEnergy - Ev) / max(gasMixtureCv,SMALL);

            if (TkLimit > 0.9 * tvTemperature)
            {
                TkLimit = 0.9 * tvTemperature;
            }
            else if(TkLimit < -9.0 * tvTemperature)
            {
                TkLimit = -9.0 * tvTemperature;
            }
            tvTemperature = MAX(tvTemperature - TkLimit, veTemperatureMinNonDimensional);

            if (fabs(TkLimit / tvTemperature) <= precisionTemperature)
            {
                break;
            }
        }
        //tvTemperature = MAX(tvTemperature, veTemperatureMinNonDimensional);

        transRotationTemperature = (pressure - density * gasConstant * electronTemperature * ceDivideMe) / (density * gasConstant * massReciprocal);
        transRotationTemperature = MAX(transRotationTemperature, SMALL);
    }

   return;



    RDouble  electronPressure;
    if (ntmodel == 2)      //! Two-Temperature model.
    {
        RDouble Tve = MAX(primitiveVariables[nl + nchem], 30.0/referenceTemperature);
        primitiveVariables[nl + nchem] = GetMixedGasVibrationEnergy(Tve, massFractions) + GetMixedGasElectronEnergy(Tve, massFractions);

        //! Compute pressure of electron.
        electronPressure = density * gasConstant * Tve * ceDivideMe;
            
        transRotationTemperature = (pressure - electronPressure) / (density * gasConstant * massReciprocal);
        vibrationTemperature = Tve;
    }
    else  //! Three-Temperature model.
    {
        RDouble Tv = primitiveVariables[nl + nchem];
        primitiveVariables[nl + nchem] = GetMixedGasVibrationEnergy(Tv, massFractions);

        RDouble Te = primitiveVariables[nl + nchem + 1];
        primitiveVariables[nl + nchem + 1] = GetMixedGasElectronEnergy(Te, massFractions);

        //! Compute pressure of electron.
        electronPressure = density * gasConstant * Te * ceDivideMe;
    
        transRotationTemperature = (pressure - electronPressure) / (density * gasConstant * massReciprocal);
        vibrationTemperature = Tv;
        electronTemperature = Te;
    }
}

RDouble MixedGas::GetTemperatureMin()
{
    RDouble refDimensionalTemperature = this->referenceParameterDimensional.GetTemperature();
    return this->thermodynamicManager->GetTemperatureMin() / refDimensionalTemperature;
}

RDouble MixedGas::GetTemperatureMax()
{
    RDouble refDimensionalTemperature = this->referenceParameterDimensional.GetTemperature();
    return this->thermodynamicManager->GetTemperatureMax() / refDimensionalTemperature;
}

string * MixedGas::GetNameOfSpecies() const 
{
    if ( !molecularProperty ) return 0;
    return molecularProperty->GetNameOfSpecies(); 
};

void MixedGas::GetElectronMolecularWeightAndMassFraction(RDouble *massFraction, RDouble &ce, RDouble &Me)
{
    RDouble Mref = referenceParameterDimensional.GetAverageMolecularWeight();
    //! Initialize.
    ce = 0.0;
    Me = 5.486e-7;      //! Unit: kg/mol.
    Me = Me / Mref;
    if (nElectronIndex >= 0)
    {
        RDouble *Ms = GetMolecularWeight();
        ce = massFraction[nElectronIndex];
        Me = Ms[nElectronIndex];
    }
}

void MixedGas::GetPressure(RDouble *primitiveVars, RDouble gama, RDouble internalEnergy, RDouble temperature, RDouble &pressure)
{
    if (nchem != 1)      //! Perfect gas.
    {
        pressure = (gama - one) * internalEnergy;
        return;
    }

    //int nOuterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    //int nChemicalFlowStep = GlobalDataBase::GetIntParaFromDB("nChemicalFlowStep");

    using namespace IDX;
    RDouble density, specificInternalEnergy, gasConstant;
    gasConstant = this->GetUniversalGasConstant();
    density = primitiveVars[IR];
    specificInternalEnergy = internalEnergy / density;

    //! Obtain the mass fractions of species.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVars[nm + i];
    }

    RDouble transRotationEnergy, transRotationTemperature, electronPressure, massReciprocal, ceDivideMe, Tk;
    massReciprocal = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFractions, ceDivideMe);
    Tk = temperature;
    if (ntmodel == 2)     //! Two-Temperature model.
    {
        RDouble vibrationElectronEnergy = primitiveVars[nl + nchem];
        RDouble vibrationElectronTemperature = this->ComputeVibrationAndElectronTemperature(massFractions, vibrationElectronEnergy, Tk);

        transRotationEnergy = specificInternalEnergy - vibrationElectronEnergy;
        transRotationTemperature = this->ComputeTranslationAndRotationTemperature(massFractions, transRotationEnergy);
        electronPressure = density * vibrationElectronTemperature * gasConstant * ceDivideMe;

        pressure = density * transRotationTemperature * gasConstant * massReciprocal + electronPressure;
    }
    else if (ntmodel == 3)     //! Three-Temperature model.
    {
        RDouble electronEnergy = primitiveVars[nl + nchem + 1];
        RDouble electronTemperature = this->ComputeElectronTemperature(massFractions, electronEnergy, Tk);

        transRotationEnergy = specificInternalEnergy - primitiveVars[nl + nchem] - electronEnergy;
        transRotationTemperature = this->ComputeTranslationAndRotationTemperature(massFractions, transRotationEnergy);
        electronPressure = density * electronTemperature * gasConstant * ceDivideMe;

        pressure = density * transRotationTemperature * gasConstant * massReciprocal + electronPressure;
    }
    else       //! The default is One-Temperature model.
    {
        transRotationTemperature = ComputeNonequilibrumTemperature(massFractions, specificInternalEnergy, Tk);
        pressure = density * gasConstant * transRotationTemperature * (massReciprocal + ceDivideMe);
    }
}

void MixedGas::GetTemperatureAndPressure(RDouble *primitiveVars, RDouble internalEnergy, RDouble *temperature, RDouble &pressure)
{
    using namespace IDX;
    RDouble density, specificInternalEnergy, gasConstant;
    gasConstant = this->GetUniversalGasConstant();
    density = primitiveVars[IR];
    specificInternalEnergy = internalEnergy / density;

    //! Obtain the mass fractions of species.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVars[nm + i];
    }

    RDouble transRotationEnergy, transRotationTemperature, electronPressure, massReciprocal, ceDivideMe, Tk;
    massReciprocal = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFractions, ceDivideMe);
    
    if (ntmodel == 2)     //! Two-Temperature model.
    {
        Tk = temperature[ITV];
        RDouble vibrationElectronEnergy = primitiveVars[nl + nchem];
        RDouble vibrationElectronTemperature = this->ComputeVibrationAndElectronTemperature(massFractions, vibrationElectronEnergy, Tk);

        transRotationEnergy = specificInternalEnergy - vibrationElectronEnergy;
        transRotationTemperature = this->ComputeTranslationAndRotationTemperature(massFractions, transRotationEnergy);
        electronPressure = density * vibrationElectronTemperature * gasConstant * ceDivideMe;

        pressure = density * transRotationTemperature * gasConstant * massReciprocal + electronPressure;

        //Save the temperature of the next time-advancing step.
        temperature[ITT] = transRotationTemperature;
        temperature[ITV] = vibrationElectronTemperature;
        temperature[ITE] = temperature[ITV];
    }
    else if (ntmodel == 3)     //! Three-Temperature model.
    {
        Tk = temperature[ITE];
        RDouble electronEnergy = primitiveVars[nl + nchem + 1];
        RDouble electronTemperature = this->ComputeElectronTemperature(massFractions, electronEnergy, Tk);

        transRotationEnergy = specificInternalEnergy - primitiveVars[nl + nchem] - electronEnergy;
        transRotationTemperature = this->ComputeTranslationAndRotationTemperature(massFractions, transRotationEnergy);
        electronPressure = density * electronTemperature * gasConstant * ceDivideMe;

        pressure = density * transRotationTemperature * gasConstant * massReciprocal + electronPressure;

        //Compute the vibrational temperature.
        RDouble vibrationEnergy = primitiveVars[nl + nchem];
        Tk = temperature[ITV];
        RDouble vibrationTemperature = this->ComputeVibrationTemperature(massFractions, vibrationEnergy, Tk);

        //Save the temperature of the next time-advancing step.
        temperature[ITT] = transRotationTemperature;
        temperature[ITV] = vibrationTemperature;
        temperature[ITE] = electronTemperature;
    }
    else       //! The default is One-Temperature model.
    {
        Tk = temperature[ITT];
        transRotationTemperature = ComputeNonequilibrumTemperature(massFractions, specificInternalEnergy, Tk);
        pressure = density * gasConstant * transRotationTemperature * (massReciprocal + ceDivideMe);

        //Save the temperature of the next time-advancing step.
        temperature[ITT] = transRotationTemperature;
        temperature[ITV] = transRotationTemperature;
        temperature[ITE] = transRotationTemperature;
    }
}

void MixedGas::NormalizePrimitive(RDouble *primitiveVariables)
{
    RDouble fSum = 0.0;
    //! Obtain the index of electron.
    if (this->nElectronIndex < 0)      //! No electron.
    {
        for (int s = 0; s < this->nLastSpeciesIndex; ++ s)
        {
            fSum += primitiveVariables[nm + s];
        }

        //! Compute mass fraction of the last species.
        primitiveVariables[nm + nLastSpeciesIndex] = MAX(0.0, 1.0 - fSum);
        fSum += primitiveVariables[nm + nLastSpeciesIndex];
    }
    else       //! The electron exists. the ions must neutralize.
    {
        int *ionType = molecularProperty->GetIonTypeOfSpecies();
        RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocal();
        primitiveVariables[nm + nElectronIndex] = 0.0;    //! Electron.
        for (int s = 0; s < this->nLastSpeciesIndex; ++ s)
        {
            primitiveVariables[nm + nElectronIndex] += ionType[s] * speciesMassReciprocal[s] * primitiveVariables[nm + s];
            fSum += primitiveVariables[nm + s];
        }
        primitiveVariables[nm + nElectronIndex] /= speciesMassReciprocal[nElectronIndex];
        primitiveVariables[nm + nLastSpeciesIndex] = MAX(1.0 - primitiveVariables[nm + nElectronIndex] - fSum, 0.0);
        fSum += primitiveVariables[nm + nElectronIndex] + primitiveVariables[nm + nLastSpeciesIndex];
    }

    //! Normalize the mass fractions.
    for (int s = 0; s < this->numberOfSpecies; ++ s)
    {
        primitiveVariables[nm + s] /= fSum;
    }
}

void MixedGas::ComputeMolecularWeightReciprocal(RDouble *primitiveVariables, RDouble &massReciprocal)
{
    RDouble *oMolecularWeight = molecularProperty->GetMolecularWeightReciprocal();

    massReciprocal = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        massReciprocal += primitiveVariables[ispecies + nm] * oMolecularWeight[ispecies];
    }
}

void MixedGas::ComputeMolecularWeightReciprocalDimensional(RDouble *prim, RDouble &omavDimensional)
{
    RDouble *oMolecularWeightDimensional = molecularProperty->GetMolecularWeightReciprocalDimensional();

    omavDimensional = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        omavDimensional += prim[ispecies + nm] * oMolecularWeightDimensional[ispecies];
    }
}

void MixedGas::ComputeOneTemperatureModelSpeciesViscosity(RDouble temperature, RDouble *viscosityOfSpecies)
{
    //! Obtain the dimensional values of each species viscosity.
    ComputeOneTemperatureModelSpeciesViscosityDimensional(temperature, speciesViscosity);

    RDouble referenceViscosity = this->referenceParameterDimensional.GetViscosity();

    //! Translate to the non-dimensional value.
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        viscosityOfSpecies[ispecies] = speciesViscosity[ispecies] / referenceViscosity;
    }
}

void MixedGas::ComputeOneTemperatureModelSpeciesViscosityDimensional(RDouble temperature, RDouble *viscosityOfSpeciesDimensional)
{
    //! To obtain the dimensional viscosity of each species according to the Blotter fitting formula.
    //! temperature is the non-dimensional temperature of the mixture gas.
    //! viscosityOfSpeciesDimensional is an array saving the dimensional viscosity of each species [N*s/m2].

    int nType = 0;     //! 0 stands for selecting the Blotter curve fits mode, 1 is for Sutherland relation.
    nType = GlobalDataBase::GetIntParaFromDB("nSutherland");
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble T = temperature * referenceTemperature;

    if (nType == 0)     //! Compute species viscosity using the Blotter curve fits for single-temperature model.
    {
        RDouble *blotterAi = CurveFits->GetBlotterCurveFittingCoefficientA();
        RDouble *blotterBi = CurveFits->GetBlotterCurveFittingCoefficientB();
        RDouble *blotterCi = CurveFits->GetBlotterCurveFittingCoefficientC();
        RDouble *blotterDi = CurveFits->GetBlotterCurveFittingCoefficientD();
        RDouble *blotterEi = CurveFits->GetBlotterCurveFittingCoefficientE();
   
        RDouble lnt1, lnt2, lnt3, lnt4, tmp;
        lnt1 = log(T);
        lnt2 = lnt1 * lnt1;
        lnt3 = lnt1 * lnt2;
        lnt4 = lnt1 * lnt3;

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            tmp = blotterAi[ispecies] * lnt4
                + blotterBi[ispecies] * lnt3
                + blotterCi[ispecies] * lnt2
                + blotterDi[ispecies] * lnt1
                + blotterEi[ispecies];

            //! The Blotter curve fits.
            viscosityOfSpeciesDimensional[ispecies] = 0.10 * exp(tmp) + SMALL;
        }
    }
    else       //! Compute species viscosity using the Sutherland relation.
    { 
        RDouble miu0, T0, C;
        miu0 = 1.7894E-5;     //! Unit: kg/(m*s).
        T0 = 288.15;          //! Unit: K.
        C  = 110.4;           //! Unit: K.
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            viscosityOfSpeciesDimensional[ispecies] = miu0 * pow(T / T0, 1.5) * (T0 + C)/(T + C);
        }
    }
}

void MixedGas::ComputeSpeciesEnthalpy(RDouble temperature, RDouble *speciesCv, RDouble *speciesE, RDouble *speciesH)
{
    ComputeSpeciesConstantVolumeSpecificHeat(temperature, speciesCv);

    ComputeSpeciesConstantPressureSpecificHeat(temperature, speciesH);

    RDouble gasConstant = this->GetCoefficientOfStateEquation();
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        speciesE[s] = speciesH[s] - gasConstant / speciesMass[s];
    }
}

void MixedGas::ComputeSpeciesConstantPressureSpecificHeatByEuckenFormula(RDouble *cpSpecies)
{
    RDouble *molecularWeightDimensional = molecularProperty->GetMolecularWeightDimensional();

    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();

    RDouble dimension = referenceVelocity * referenceVelocity / referenceTemperature;

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        workSpecies[ispecies] = 1.0;
    }

    workSpecies[0] = 0.0;
    workSpecies[3] = 0.0;
    workSpecies[5] = 0.0;

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        cpSpecies[ispecies] = (2.5 * 1.5 + workSpecies[ispecies])
                            * (rjmk / molecularWeightDimensional[ispecies]) / dimension;
    }

    return;
}

void MixedGas::ComputeEnthalpyByPrimitive(RDouble *primitiveVariables, RDouble gama, RDouble &enthalpy, RDouble* temperatures)
{
    using namespace IDX;

    RDouble &density  = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];
    int nIdealState = GlobalDataBase::GetIntParaFromDB("nIdealState");
    if (nchem == 0 || nIdealState ==1)
    //if (nchem == 0)
    {
        enthalpy = (gama / (gama - one)) * (pressure / density);
    }
    else
    {
        enthalpy = this->ComputeMixedGasEnthalpy(primitiveVariables, temperatures);
    }
}

void MixedGas::ComputeInternalEnergy(RDouble *primitiveVariables, RDouble gama, RDouble vibrationTemperature, RDouble electronTemperature, RDouble &totalEnergy)
{
    using namespace IDX;
    RDouble temperature = 0.0, squareVelocity = 0.0, massReciprocal = 0.0, staticEnthalpy = 0.0;

    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];
    squareVelocity = um * um + vm * vm + wm * wm;
    int nIdealState = GlobalDataBase::GetIntParaFromDB("nIdealState");
    if (nchem == 0 || nIdealState ==1)//! Perfect gas.
    //if (nchem == 0)
    {
        totalEnergy = (pressure / density) / (gama - one) + half * squareVelocity;
    }
    else
    {
        if (ntmodel == 1)     //! Single-Temperature model.
        {
            ComputeMolecularWeightReciprocal(primitiveVariables, massReciprocal);
            temperature = pressure / (coefficientOfStateEquation * density * massReciprocal);
            ComputeSpeciesEnthalpy(temperature, speciesEnthalpy);
            ComputeMixtureByPrimitive(primitiveVariables, speciesEnthalpy, staticEnthalpy);
            totalEnergy = staticEnthalpy - pressure / density + half * squareVelocity;
        }
        else      //! Multi-temperature model.
        {
            for (int i = 0; i < this->numberOfSpecies; ++ i)
            {
                massFractions[i] = primitiveVariables[nm + i];
            }

            RDouble gasConstant = GetUniversalGasConstant();
            //! To obtain the static enthalpy of the mixed gas.
            RDouble electronPressure = 0.0, Ttr = 0.0, ceDivideMe = 0.0;
            massReciprocal = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFractions, ceDivideMe);

            if (ntmodel == 2)      //! Two-Temperature Model.
            {
                electronPressure = density * gasConstant * vibrationTemperature * ceDivideMe;
                //! To obtain the translation-rotation temperature where the pressure is known.
                Ttr = (pressure - electronPressure) / (density * gasConstant * massReciprocal);
                //! Compute the static enthalpy.
                staticEnthalpy = ComputeMixedGasEnthalpy(Ttr, vibrationTemperature, vibrationTemperature, massFractions);
            }
            else if (ntmodel == 3)  //! Three-Temperature Model.
            {
                electronPressure = density * gasConstant * electronTemperature * ceDivideMe;
                //! To obtain the translation-rotation temperature where the pressure is known.
                Ttr = (pressure - electronPressure) / (density * gasConstant * massReciprocal);
                
                //! Compute the static enthalpy.
                staticEnthalpy = ComputeMixedGasEnthalpy(Ttr, vibrationTemperature, electronTemperature, massFractions);
            }

            //! Compute the total internal energy.
            totalEnergy = staticEnthalpy - pressure / density + half * squareVelocity;
        }
    }
}

void MixedGas::ComputeDensityDiffusionAndKCP(RDouble temperature, RDouble *primface, RDouble mul, RDouble mut, RDouble *rho_ds_face, RDouble *hintSpeciesOfFace, RDouble &kcp)
{
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        massFractions[ispecies] = primface[nm + ispecies];
    }

    ComputeMoleFractionByMassFraction(massFractions, moleFractions);

    RDouble *oscl = GetLaminarSchmidtNumberReciprocal();
    RDouble *osct = GetTurbulentSchmidtNumberReciprocal();

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        RDouble tmp = 1.0 - moleFractions[ispecies];
        if (tmp <= 0.0)
        {
            rho_ds_face[ispecies] = mul * oscl[ispecies] + mut * osct[ispecies];
        }
        else
        {
            rho_ds_face[ispecies] = (mul * oscl[ispecies] + mut * osct[ispecies]) * (1.0 - massFractions[ispecies]) / tmp;
        }
    }


    ComputeSpeciesEnthalpy(temperature, speciesEnthalpy);

    RDouble cp = -1.0;
    ComputeConstantPressureSpecificHeatByMassFraction(massFractions, temperature, cp);
    kcp = (mul * oprl + mut * oprt) * cp;

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        hintSpeciesOfFace[ispecies] = speciesEnthalpy[ispecies];
    }
}

void MixedGas::ComputeDHAndTotalEnthalpy(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &deltaEnthalpy, RDouble &totalEnthalpy)
{
    using namespace IDX;

    RDouble squareVelocity, ae, af;
    RDouble enthalpy;
    //! Obtain the primitive variables.
    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];

    squareVelocity = um * um + vm * vm + wm * wm;
    ae = gama - one;
    af = half * ae * squareVelocity;
    deltaEnthalpy = -ae * (um * deltaQ[IRU] + vm * deltaQ[IRV] + wm * deltaQ[IRW] - deltaQ[IRE]);

    int nIdealState = GlobalDataBase::GetIntParaFromDB("nIdealState");
    if (nchem == 0 || nIdealState ==1)
    //if (nchem == 0)
    {
        //! h = Cp*T = (gamma*Rg/gamma-1)*T = (gamma/gamma-1)*(p/rho).
        enthalpy = (gama / (gama - one)) * (pressure / density);
    }
    else
    {
        //! The obtained value is the reciprocal value of molecular weight.
        RDouble *oMolecularWeight = molecularProperty->GetMolecularWeightReciprocal();

        RDouble temperature, massReciprocal, as;
        ComputeMolecularWeightReciprocal(primitiveVariables, massReciprocal);
        temperature = pressure / (coefficientOfStateEquation * density * massReciprocal);
        ComputeSpeciesEnthalpy(temperature, speciesEnthalpy);
        ComputeMixtureByPrimitive(primitiveVariables, speciesEnthalpy, enthalpy);

        for (int ispecies = 0; ispecies < numberOfSpecies - 1; ++ ispecies)
        {
            as = -ae * (speciesEnthalpy[ispecies] - speciesEnthalpy[numberOfSpecies - 1])
               + coefficientOfStateEquation * gama * temperature * (oMolecularWeight[ispecies] - oMolecularWeight[numberOfSpecies - 1]);
            deltaEnthalpy += as * deltaQ[nm + ispecies];
        }
        af += coefficientOfStateEquation * gama * temperature * oMolecularWeight[numberOfSpecies - 1] - ae * speciesEnthalpy[numberOfSpecies - 1];
    }

    deltaEnthalpy += af * deltaQ[IR];
    // ht = h+v*v/2
    totalEnthalpy = enthalpy + half * squareVelocity;
}

void MixedGas::ComputeTotalEnthalpyAndDH(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &totalEnthalpy, RDouble &deltaEnthalpy)
{
    //! This function is used for perfect gas and single temperature model.
    //! To compute the total enthalpy and the variable dh in the function MXDQ_STD(). dh=b2, b2 denotes the coefficient \n
    //! of the vector M*dQ which can be referred to the forumla (A.7) and (A.8) in the appendix A of the PHengLEI Theory manual.
    using namespace IDX;

    //! Obtain the primitive variables.
    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];

    //! c1 is product coefficient of the variable dq[IR]£¬c2 = R*T/Mns - (gama-1)*ens£¬gama1 = gama-1.
    RDouble squareVelocity, c1, c2, gama1;
    squareVelocity = um * um + vm * vm + wm * wm;
    gama1 = gama - 1.0;
    c1 = half * gama1 * squareVelocity;
    deltaEnthalpy = -gama1 * (um * deltaQ[IRU] + vm * deltaQ[IRV] + wm * deltaQ[IRW] - deltaQ[IRE]);

    //! enthalpy is static enthalpy.
    RDouble enthalpy;
    int nIdealState = GlobalDataBase::GetIntParaFromDB("nIdealState");
    if (nchem == 0 || nIdealState == 1)     //! Compute the static enthalpy.
    //if (nchem == 0)
    {
        // h = cp*T = [gama/(gama-1)]*(R/M)*T = [gama/(gama-1)]*(p/rho).
        enthalpy = (gama / gama1)* (pressure / density);
    }
    else      //! Compute the static enthalpy of the multi-species model.
    {
        //! Obtain the non-dimensional molecular weight of each species saving in the array Ms.
        RDouble *speciesMass = molecularProperty->GetMolecularWeight();
        RDouble gasConstant, temperature, mass;
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            massFractions[i] = primitiveVariables[nm + i];
        }

        gasConstant = GetUniversalGasConstant();
        mass = GetMixedGasMolecularWeight(massFractions);
        //! Compute temperature using the state equation.
        temperature = (pressure / density) * (mass / gasConstant);

        //! Obtain constant volume specific heat of each species.
        ComputeSpeciesConstantVolumeSpecificHeat(temperature, this->speciesCv);

        // c2 = R*T/Mns - (gama-1)*ens
        c2 = gasConstant * temperature / speciesMass[numberOfSpecies - 1] - gama1 * speciesCv[numberOfSpecies - 1] * temperature;
        c1 += c2;     //! Added to the product coefficient of the variable dq[IR].

        //! fSum saves the summation of each term.
        RDouble fSum = 0.0;
        for (int i = 0; i < numberOfSpecies - 1; ++ i)
        {
            fSum += deltaQ[nm + i];
        }
        deltaEnthalpy += -c2 * fSum;

        fSum = 0.0;
        for (int i = 0; i < numberOfSpecies - 1; ++ i)
        {
            fSum += deltaQ[nm + i] / speciesMass[i];
        }
        deltaEnthalpy += gasConstant * temperature * fSum;

        fSum = 0.0;
        for (int i = 0; i < numberOfSpecies - 1; ++ i)
        {
            fSum += deltaQ[nm + i] * speciesCv[i] * temperature;
        }
        deltaEnthalpy += -gama1 * fSum;

        this->ComputeSpeciesEnthalpy(temperature, speciesEnthalpy);
        enthalpy = 0.0;     //! Compute the static enthalpy.
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            enthalpy += massFractions[i] * speciesEnthalpy[i];
        }

    }

    deltaEnthalpy += c1 * deltaQ[IR];
    // h0 = h+v*v/2, total enthalpy equals to static enthalpy plus half of the squared velocity.
    totalEnthalpy = enthalpy + half * squareVelocity;
}

void MixedGas::ComputeViscosityByPrimitive(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &electronTemperature, RDouble &suthTemperature, RDouble &viscosity)
{
    using namespace IDX;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        massFractions[ispecies] = primitiveVariables[nm + ispecies];
    }

    ComputeViscosity(primitiveVariables[IR], massFractions, transRotationTemperature, electronTemperature, suthTemperature, viscosity);
}

void MixedGas::ComputeTransportCoefficients(RDouble *primitiveVariables, RDouble *temperature, RDouble visTurbulence, RDouble &viscosity, RDouble *heatConductivity, RDouble *rhoDs)
{
    using namespace IDX;
    //RDouble rho = primitiveVariables[IR];
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        massFractions[s] = primitiveVariables[nm + s];
    }
    // Compute mole fractions using the mass fraction.
    ComputeMoleFractionByMassFraction(massFractions, moleFractions);

    RDouble *speciesMass = molecularProperty->GetMolecularWeight();

    if (ntmodel > 1) //Multi-temperature model.
    {
        ComputeSpeciesViscosityWithCurveFitMethod(temperature[ITT], temperature[ITE], speciesViscosity);

        GetPartitionFunctionValues(speciesViscosity, speciesMass, moleFractions, speciesWeight);

        // Compute heat conductivity.
        ComputeSpeciesTranslationAndRotationHeatConductivity(speciesViscosity, speciesConductivity);
        ComputeMixtureByWilkeFormula(moleFractions, speciesConductivity, speciesWeight, heatConductivity[ITT]);

        ComputeVibrationSpecificHeatAtConstantVolume(temperature[ITV], speciesCv);
        ComputeSpeciesVibrationHeatConductivity(speciesViscosity, speciesCv, speciesConductivity);
        ComputeMixtureByWilkeFormula(moleFractions, speciesConductivity, speciesWeight, heatConductivity[ITV]);

        ComputeElectronSpecificHeatAtConstantVolume(temperature[ITE], speciesCv);
        ComputeSpeciesElectronHeatConductivity(speciesViscosity, speciesCv, speciesConductivity);
        ComputeMixtureByWilkeFormula(moleFractions, speciesConductivity, speciesWeight, heatConductivity[ITE]);
    }
    else//One-temperature model.
    {
        ComputeOneTemperatureModelSpeciesViscosity(temperature[ITT], speciesViscosity);
        GetPartitionFunctionValues(speciesViscosity, speciesMass, moleFractions, speciesWeight);

        ComputeSpeciesConstantPressureSpecificHeat(temperature[ITT], speciesCp);
        ComputeSpeciesHeatConductivityByEuckenFormula(speciesViscosity, speciesCp, speciesConductivity);

        // Compute heat conductivity.
        ComputeMixtureByWilkeFormula(moleFractions, speciesConductivity, speciesWeight, heatConductivity[ITT]);
    }

    // Compute viscosity.
    ComputeMixtureByWilkeFormula(moleFractions, speciesViscosity, speciesWeight, viscosity);

    // Compute mass diffusion coefficients.
    ComputeSpeciesMassDiffusionCoefficient(massFractions, viscosity, visTurbulence, rhoDs);

}

void MixedGas::ComputeTransportCoefficientsR(RDouble *massFractions, RDouble *temperature, RDouble *speciesCps, RDouble *speciesCvvs,RDouble *speciesCves, RDouble visTurbulence, RDouble &viscosity, RDouble *heatConductivity, RDouble *rhoDs)
{
    using namespace IDX;

    RDouble *oMolecularWeight = molecularProperty->GetMolecularWeightReciprocal();     //! The reciprocal value of mass fraction.
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceViscosity = this->referenceParameterDimensional.GetViscosity();
    int *ionType = this->molecularProperty->GetIonTypeOfSpecies();
    RDouble *oSchmidtLaminar = GetLaminarSchmidtNumberReciprocal();
    RDouble *oSchmidtTurbulence = GetTurbulentSchmidtNumberReciprocal();
    RDouble gasConstant = this->GetUniversalGasConstant(); 

    RDouble *blotterAi = CurveFits->GetBlotterCurveFittingCoefficientA();
    RDouble *blotterBi = CurveFits->GetBlotterCurveFittingCoefficientB();
    RDouble *blotterCi = CurveFits->GetBlotterCurveFittingCoefficientC();
    RDouble *blotterDi = CurveFits->GetBlotterCurveFittingCoefficientD();
    RDouble *blotterEi = CurveFits->GetBlotterCurveFittingCoefficientE();

    RDouble Ttr = temperature[mTT] * referenceTemperature;

    RDouble xiSum, oxiSum;
    xiSum = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        moleFractions[ispecies] = massFractions[ispecies] * oMolecularWeight[ispecies];
        xiSum += moleFractions[ispecies];
    }

    oxiSum = 1.0 / xiSum;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
       moleFractions[ispecies] *= oxiSum;
    }

    RDouble lnT = log(static_cast<RealData>(Ttr));
    for (int ispecies = 0; ispecies <= nLastSpeciesIndex; ++ ispecies)
    {
        oxiSum = min(maxViscous,(((blotterAi[ispecies] * lnT + blotterBi[ispecies]) * lnT + blotterCi[ispecies]) * lnT + blotterDi[ispecies]) * lnT + blotterEi[ispecies]);
        speciesViscosity[ispecies] = 0.10 * exp(static_cast<RealData>(oxiSum)) / referenceViscosity + SMALL;
    }

    if(nElectronIndex >= 0)
    {
        RDouble Te  = temperature[mTE] * referenceTemperature;
        RDouble lnTe = log(static_cast<RealData>(Te));
        oxiSum = min(maxViscous,(((blotterAi[nElectronIndex] * lnTe + blotterBi[nElectronIndex]) * lnTe + blotterCi[nElectronIndex]) * lnTe + blotterDi[nElectronIndex]) * lnTe + blotterEi[nElectronIndex]);
        speciesViscosity[nElectronIndex] = 0.10 * exp( static_cast<RealData>(oxiSum)) / referenceViscosity + SMALL;
    }
    
    RDouble tmp0, tmp1, phiTmp;
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        phiTmp = 0.0;
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            tmp0 = speciesMass[i] / speciesMass[s];
            tmp1 = 1.0 + sqrt(static_cast<RealData>(speciesViscosity[s] / speciesViscosity[i])) * sqrt((sqrt(static_cast<RealData>(tmp0))));
            phiTmp += moleFractions[i] * tmp1 * tmp1 / sqrt(static_cast<RealData>(8.0 + 8.0 / tmp0));
        }
        speciesWeight[s] = phiTmp;
    }

    //! Compute the viscosity using Eucken formula.
    viscosity = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        speciesConductivity[ispecies] = moleFractions[ispecies] * speciesViscosity[ispecies] / speciesWeight[ispecies];
        viscosity += speciesConductivity[ispecies];
    }

    //! Compute mass diffusion coefficient for Species
    tmp1 = 0;
    phiTmp = 0;
    for(int i = 0; i <= this->nLastSpeciesIndex; ++ i)
    {
        if (1.0 - moleFractions[i] <= 1.0e-20)
        {
            rhoDs[i] = 0.0;
        }
        else
        {
            rhoDs[i] = (1.0 - massFractions[i]) / (1.0 - moleFractions[i]) * (viscosity * oSchmidtLaminar[i] + visTurbulence * oSchmidtTurbulence[i]);
        }

        if (ionType[i] != 0)
        {
            tmp0  = massFractions[i] * oMolecularWeight[i];
            tmp1 += rhoDs[i] * tmp0;
            phiTmp += tmp0;
        }
    }

    //! Compute the mass diffusion coefficient of electron.
    if (nElectronIndex >= 0)
    {
        rhoDs[nElectronIndex] = tmp1 / MAX(phiTmp, 1.0e-20);
    }

    //! Compute the heat conductivity using Eucken formula.
    heatConductivity[ITT] = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        heatConductivity[ITT] += speciesConductivity[ispecies] * ( speciesCps[ispecies] + 1.25 * gasConstant * oMolecularWeight[ispecies]);
    }

    heatConductivity[ITV] = 0.0;
    heatConductivity[ITE] = 0.0;
    if(ntmodel == 2)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            heatConductivity[ITV] += speciesConductivity[ispecies] * ( speciesCvvs[ispecies] + speciesCves[ispecies] );
        }
    }
    else if (ntmodel == 3)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            heatConductivity[ITV] += speciesConductivity[ispecies] * speciesCvvs[ispecies];
            heatConductivity[ITE] += speciesConductivity[ispecies] * speciesCves[ispecies];
        }
    }
}


void MixedGas::ComputeViscosity(RDouble density, RDouble *massFractionOfSpecies, RDouble &transRotationTemperature, RDouble &electronTemperature, RDouble &suthTemperature, RDouble &viscosity)
{
    if (nchem == 0)     //! Perfect gas.
    {
        viscosity = transRotationTemperature * sqrt(transRotationTemperature) * (1.0 + suthTemperature) / (transRotationTemperature + suthTemperature);
    }
    else      //! Chemical non-equilibrium flow.
    {
        RDouble *speciesMass = molecularProperty->GetMolecularWeight();
        ComputeMoleFractionByMassFraction(massFractionOfSpecies, moleFractions);

        if (ntmodel > 1)     //! Multi-temperature model.
        {
            RDouble electronPressure = 0.0;
            int nType = 0;    //! 0 for Blottner curve fit method, 1 for Lennard-Jones method.
            if (nType == 1)
            {
                double electronMassFraction = 0.0;
                if (nElectronIndex >= 0)
                {
                    electronMassFraction = massFractionOfSpecies[nElectronIndex];
                }
                RDouble gasConstant = GetUniversalGasConstant();
                electronPressure = (density * electronTemperature * gasConstant) * (electronMassFraction / electronMass);
                ComputeSpeciesViscosityWithLennardJonesMethod(transRotationTemperature, electronPressure, speciesViscosity);
            }
            else
            {
                ComputeSpeciesViscosityWithCurveFitMethod(transRotationTemperature, electronTemperature, speciesViscosity);
            }
        }
        else     //! One-temperature model.
        {
            ComputeOneTemperatureModelSpeciesViscosity(transRotationTemperature, speciesViscosity);
        }

        GetPartitionFunctionValues(speciesViscosity, speciesMass, moleFractions, speciesWeight);

        viscosity = GetMixedGasViscosityWithWilkeFormula(moleFractions, speciesViscosity, speciesWeight);
    }
}

void MixedGas::ComputeSpeciesHeatConductivityByEuckenFormula(RDouble *viscosity, RDouble *cps, RDouble *conductivity)
{
    //! Obtain the universal gas constant.
    RDouble gasConstant = this->GetUniversalGasConstant(); 
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();

    //! Compute the heat conductivity using Eucken formula.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        conductivity[i] = viscosity[i] * cps[i];
        conductivity[i] += 1.25 * viscosity[i] * gasConstant / speciesMass[i];
    }
}

RDouble MixedGas::ComputeMixtureGasHeatConductivityByWassilewaFormula(RDouble *primitiveVariables, RDouble temperature)
{
    //int nSpecies = this->GetNumberOfSpecies();

    //! Allocate the memory spaces.
    //RDouble *massFraction = new RDouble[nSpecies];
    //RDouble *moleFraction = new RDouble[nSpecies];// To store the mole fraction of each species.
    //RDouble *speciesCp = new RDouble[nSpecies];   // To store the species' specific heat at constant pressure.
    //RDouble *phi = new RDouble[nSpecies];         // To store the partition function of each species.
    //RDouble *viscosity = new RDouble[nSpecies];   // To store the viscosity of each species.
    //RDouble *conductivity = new RDouble[nSpecies];// To store the heat conductivity of each species.

    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVariables[nm + i];
    }

    //! Obtain the mole fraction of each species.
    this->ComputeMoleFractionByMassFraction(massFractions, moleFractions);

    //! Obtain the molecular weight of each species.
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();

    //! Obtain the viscosity of each species.
    this->ComputeOneTemperatureModelSpeciesViscosity(temperature, speciesViscosity);

    //! Obtain the species specific heat ratio at constant pressure.
    this->ComputeSpeciesConstantPressureSpecificHeat(temperature, speciesCp);

    //! Compute the heat conductivity of each species.
    this->ComputeSpeciesHeatConductivityByEuckenFormula(speciesViscosity, speciesCp, speciesConductivity);

    //! Compute species partition function values.
    this->GetPartitionFunctionValues(speciesViscosity, speciesMass, moleFractions, this->workSpecies);

    RDouble lambda = 0.0;
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        lambda += speciesConductivity[i] * moleFractions[i] / workSpecies[i];
    }

    //! Remove the allocated memory space.
    //delete []massFraction;
    //delete []moleFraction;
    //delete []speciesCp;
    //delete []phi;
    //delete []viscosity;
    //delete []conductivity;
    // Return the value.
    return lambda;
}

void MixedGas::GetPartitionFunctionValues(RDouble *target, RDouble *speciesMass, RDouble *moleFraction, RDouble *phi)
{
    //! To obtain the partition function values based on the targets, e.g conductivity or viscosity.
    RDouble tmp = 0.0;
    RDouble tmp0, tmp1, phiTmp;
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        phiTmp = 0.0;
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            tmp0 = speciesMass[i] / speciesMass[s];
            tmp1 = 1.0 + sqrt(target[s] / target[i]) * sqrt(static_cast<float>(sqrt(static_cast<float>(tmp0))));
            tmp = moleFraction[i] * tmp1 * tmp1;
            phiTmp += tmp / sqrt(static_cast<float>(8.0 + 8.0 / tmp0));
        }
        phi[s] = phiTmp;
    }
}

//void RealGas::ChemicalSource(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, int nCellNumber)
//{
//    using namespace IDX;
//    RDouble *speciesMass  = molecularProperty->GetMolecularWeightDimensional();
//    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
//    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
//    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
//    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
//
//    //! The dimensional value of control temperature of chemical reaction.
//    //RDouble dimensionalDensity, forwardReactionTemperature, backwardReactionTemperature;
//    //RDouble forwardReactionItem, backwardReactionItem, collisionItem;
//    //RDouble densityPerMass, fSum;
//
//    //! The variables is used for dimensionless processing. referenceSpeciesEquation denotes the reference values of source items \n
//    //! in the species transport equations, while referenceEnergyEquation is the reference values of source items for the energy equations.
//    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
//
//    RDouble *forwardRateCoefA = reactionRate->GetForwardReactionRateCoefficientA();
//    RDouble *forwardRateCoefB = reactionRate->GetForwardReactionRateCoefficientB();
//    RDouble *forwardRateCoefC = reactionRate->GetForwardReactionRateCoefficientC();
//
//    RDouble *backwardRateCoefA = reactionRate->GetBackwardReactionRateCoefficientA();
//    RDouble *backwardRateCoefB = reactionRate->GetBackwardReactionRateCoefficientB();
//    RDouble *backwardRateCoefC = reactionRate->GetBackwardReactionRateCoefficientC();
//
//    int *reactionType = stoichiometricEquation->GetChemicalReactionType();
//    bool *isCollisionFlag = stoichiometricEquation->GetCollisionReactionFlag();
//
//    int **forwardStoichNumber = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
//    int **backwardStoichNumber = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
//    RDouble **thirdBodyStoichNumber = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();
//
//    RDouble *formationRates = new RDouble [numberOfReaction];     //! Formation rates of chemical reactions.
//    RDouble *forwardRates   = new RDouble [numberOfReaction];     //! Equilibrium constant of forward reaction.
//    RDouble *backwardRates  = new RDouble [numberOfReaction];     //! Equilibrium constant of backward reaction.
//
//    RDouble trsTemperature, vibTemperature, Tf[6], Tb[6];
//    for (int iCell = 0; iCell < nCellNumber; ++ iCell)
//    {
//        //! The control temperature of front reaction.
//        //RDouble forwardReactionTemperature = MAX(temperatures[ITT][iCell] * referenceTemperature, SMALL);
//
//        ////! The control temperature of back reaction. it is temperature of translation and rotation.
//        //RDouble backwardReactionTemperature = forwardReactionTemperature;
//        //if (ntmodel > 1)
//        //{
//        //    //! Referred to the Park and Sharma formula.
//        //    forwardReactionTemperature = pow(forwardReactionTemperature, 0.6) * pow(temperatures[ITV][iCell] * referenceTemperature, 0.4);
//        //}
//
//        //Obtain the dimensional temperature.
//        trsTemperature = MAX(temperatures[ITT][iCell] * referenceTemperature, SMALL);
//        vibTemperature = trsTemperature;
//        if (ntmodel > 1)
//        {
//            vibTemperature = MAX(temperatures[ITV][iCell] * referenceTemperature, SMALL);
//        }
//        //the control temperature of dissociation reaction, e.g. AB = A + B.
//        Tf[0] = pow(trsTemperature, parkVDPower) * pow(vibTemperature, 1.0 - parkVDPower);
//        Tb[0] = trsTemperature;
//        //the control temperature of electron collision dissociation reaction, e.g. AB + e- = A + B + e-.
//        Tf[1] = vibTemperature;
//        Tb[1] = sqrt(trsTemperature * vibTemperature);
//        //the control temperature of neutral exchange reaction, e.g. AB + C = A + BC.
//        Tf[2] = trsTemperature;
//        Tb[2] = trsTemperature;
//        //the control temperature of charge exchange reaction, e.g. AB+ + C = A+ + BC.
//        Tf[3] = trsTemperature;
//        Tb[3] = trsTemperature;
//        //the control temperature of recombination ionization reaction, e.g. A + B = AB+ + e-.
//        Tf[4] = trsTemperature;
//        Tb[4] = Tb[1];
//        //the control temperature of electron collision ionization reaction, e.g. A + e- = A+ + 2e-.
//        Tf[5] = vibTemperature;
//        Tb[5] = vibTemperature;
//
//        //! The control temperature of front reaction.
//        RDouble forwardReactionTemperature = Tf[0];
//        //! The control temperature of back reaction. it is temperature of translation and rotation.
//        RDouble backwardReactionTemperature = Tb[0];
//        RDouble dimensionalDensity = primitiveVariables[IR][iCell] * referenceDensity;
//        for (int iReaction = 0; iReaction < numberOfReaction; ++ iReaction)
//        {
//            RDouble forwardReactionItem = 1.0;
//            RDouble backwardReactionItem = 1.0;
//            RDouble collisionItem = 0.0;
//            int nType = reactionType[iReaction]; //Get the type of chemical reaction.
//            forwardReactionTemperature = Tf[nType];
//            backwardReactionTemperature = Tb[nType];
//
//            for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
//            {
//                RDouble densityPerMass = dimensionalDensity * primitiveVariables[nm + iSpecies][iCell] * speciesMassReciprocal[iSpecies];
//
//                for (int i = 0; i < forwardStoichNumber[iReaction][iSpecies]; ++ i)
//                {
//                    forwardReactionItem *= densityPerMass;
//                }
//                for (int i = 0; i < backwardStoichNumber[iReaction][iSpecies]; ++ i)
//                {
//                    backwardReactionItem *= densityPerMass;
//                }
//                //forwardReactionItem *= pow(densityPerMass, forwardStoichNumber[iReaction][iSpecies]);
//                //backwardReactionItem *= pow(densityPerMass, backwardStoichNumber[iReaction][iSpecies]);
//                //! The collision efficiency of species ispecies in the ireaction reaction.
//                collisionItem = collisionItem + densityPerMass * thirdBodyStoichNumber[iReaction][iSpecies]; 
//            }
//
//            if (!isCollisionFlag[iReaction])
//            {
//                collisionItem = 1.0;
//            }
//
//            //! Compute equilibrium constant of the front reaction.
//            forwardRates[iReaction] = forwardRateCoefA[iReaction] * pow(forwardReactionTemperature, forwardRateCoefB[iReaction]) * exp(-forwardRateCoefC[iReaction] / forwardReactionTemperature);
//            //! Compute equilibrium constant of the reverse reaction.
//            backwardRates[iReaction] = backwardRateCoefA[iReaction] * pow(backwardReactionTemperature, backwardRateCoefB[iReaction]) * exp(-backwardRateCoefC[iReaction] / backwardReactionTemperature);
//            //! Formation rate of the ireaction reaction.
//            formationRates[iReaction] = (forwardRates[iReaction] * forwardReactionItem - backwardRates[iReaction] * backwardReactionItem) * collisionItem;
//        }
//
//        for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
//        {
//            RDouble fSum = 0.0;
//            for (int iReaction = 0; iReaction < numberOfReaction; ++ iReaction)
//            {
//                fSum += formationRates[iReaction] * (backwardStoichNumber[iReaction][iSpecies] - forwardStoichNumber[iReaction][iSpecies]);
//            }
//            //! Compute the source terms of the specified species.
//            chemicalSourceTerms[iSpecies][iCell] = speciesMass[iSpecies] * fSum * cellVolume[iCell] * referenceSpeciesEquation;
//        }
//    }
//    delete []formationRates;
//    delete []forwardRates;
//    delete []backwardRates;
//}

void MixedGas::ChemicalSource(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, int iCell)
{
    using namespace IDX;
    RDouble *speciesMass  = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();

    //! The dimensional value of control temperature of chemical reaction.
    //RDouble dimensionalDensity, forwardReactionTemperature, backwardReactionTemperature;
    //RDouble forwardReactionItem, backwardReactionItem, collisionItem;
    //RDouble densityPerMass, fSum;

    //! The variables is used for dimensionless processing. referenceSpeciesEquation denotes the reference values of source items \n
    //! in the species transport equations, while referenceEnergyEquation is the reference values of source items for the energy equations.
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);

    RDouble *forwardRateCoefA = reactionRate->GetForwardReactionRateCoefficientA();
    RDouble *forwardRateCoefB = reactionRate->GetForwardReactionRateCoefficientB();
    RDouble *forwardRateCoefC = reactionRate->GetForwardReactionRateCoefficientC();

    RDouble *backwardRateCoefA = reactionRate->GetBackwardReactionRateCoefficientA();
    RDouble *backwardRateCoefB = reactionRate->GetBackwardReactionRateCoefficientB();
    RDouble *backwardRateCoefC = reactionRate->GetBackwardReactionRateCoefficientC();

    int *reactionType = stoichiometricEquation->GetChemicalReactionType();
    bool *isCollisionFlag = stoichiometricEquation->GetCollisionReactionFlag();

    int **forwardStoichNumber = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
    int **backwardStoichNumber = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
    RDouble **thirdBodyStoichNumber = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();

    RDouble *formationRates = new RDouble [numberOfReaction];     //! Formation rates of chemical reactions.
    RDouble *forwardRates   = new RDouble [numberOfReaction];     //! Equilibrium constant of forward reaction.
    RDouble *backwardRates  = new RDouble [numberOfReaction];     //! Equilibrium constant of backward reaction.

    RDouble trsTemperature, vibTemperature, Tf[6], Tb[6];
    trsTemperature = MAX(temperatures[ITT][iCell] * referenceTemperature, SMALL);
    vibTemperature = trsTemperature;
    if (ntmodel > 1)
    {
        vibTemperature = MAX(temperatures[ITV][iCell] * referenceTemperature, SMALL);
    }
    //the control temperature of dissociation reaction, e.g. AB = A + B.
    Tf[0] = pow(trsTemperature, parkVDPower) * pow(vibTemperature, 1.0 - parkVDPower);
    Tb[0] = trsTemperature;
    //the control temperature of electron collision dissociation reaction, e.g. AB + e- = A + B + e-.
    Tf[1] = vibTemperature;
    Tb[1] = sqrt(trsTemperature * vibTemperature);
    //the control temperature of neutral exchange reaction, e.g. AB + C = A + BC.
    Tf[2] = trsTemperature;
    Tb[2] = trsTemperature;
    //the control temperature of charge exchange reaction, e.g. AB+ + C = A+ + BC.
    Tf[3] = trsTemperature;
    Tb[3] = trsTemperature;
    //the control temperature of recombination ionization reaction, e.g. A + B = AB+ + e-.
    Tf[4] = trsTemperature;
    Tb[4] = Tb[1];
    //the control temperature of electron collision ionization reaction, e.g. A + e- = A+ + 2e-.
    Tf[5] = vibTemperature;
    Tb[5] = vibTemperature;

    //! The control temperature of front reaction.
    RDouble forwardReactionTemperature = Tf[0];
    //! The control temperature of back reaction. it is temperature of translation and rotation.
    RDouble backwardReactionTemperature = Tb[0];
    RDouble dimensionalDensity = primitiveVariables[IR][iCell] * referenceDensity;
    for (int iReaction = 0; iReaction < numberOfReaction; ++ iReaction)
    {
        RDouble forwardReactionItem = 1.0;
        RDouble backwardReactionItem = 1.0;
        RDouble collisionItem = 0.0;
        int nType = reactionType[iReaction]; //Get the type of chemical reaction.
        forwardReactionTemperature = Tf[nType];
        backwardReactionTemperature = Tb[nType];

        for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
        {
            RDouble densityPerMass = dimensionalDensity * primitiveVariables[nm + iSpecies][iCell] * speciesMassReciprocal[iSpecies];

            for (int i = 0; i < forwardStoichNumber[iReaction][iSpecies]; ++ i)
            {
                forwardReactionItem *= densityPerMass;
            }
            for (int i = 0; i < backwardStoichNumber[iReaction][iSpecies]; ++ i)
            {
                backwardReactionItem *= densityPerMass;
            }
            //forwardReactionItem *= pow(densityPerMass, forwardStoichNumber[iReaction][iSpecies]);
            //backwardReactionItem *= pow(densityPerMass, backwardStoichNumber[iReaction][iSpecies]);
            //! The collision efficiency of species ispecies in the ireaction reaction.
            collisionItem = collisionItem + densityPerMass * thirdBodyStoichNumber[iReaction][iSpecies]; 
        }

        if (!isCollisionFlag[iReaction])
        {
            collisionItem = 1.0;
        }

        //! Compute equilibrium constant of the front reaction.
        forwardRates[iReaction] = forwardRateCoefA[iReaction] * pow(forwardReactionTemperature, forwardRateCoefB[iReaction]) * exp(-forwardRateCoefC[iReaction] / forwardReactionTemperature);
        //! Compute equilibrium constant of the reverse reaction.
        backwardRates[iReaction] = backwardRateCoefA[iReaction] * pow(backwardReactionTemperature, backwardRateCoefB[iReaction]) * exp(-backwardRateCoefC[iReaction] / backwardReactionTemperature);
        //! Formation rate of the ireaction reaction.
        formationRates[iReaction] = (forwardRates[iReaction] * forwardReactionItem - backwardRates[iReaction] * backwardReactionItem) * collisionItem;
    }

    for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
    {
        RDouble fSum = 0.0;
        for (int iReaction = 0; iReaction < numberOfReaction; ++ iReaction)
        {
            fSum += formationRates[iReaction] * (backwardStoichNumber[iReaction][iSpecies] - forwardStoichNumber[iReaction][iSpecies]);
        }
        //! Compute the source terms of the specified species.
        chemicalSourceTerms[iSpecies][iCell] = speciesMass[iSpecies] * fSum * cellVolume[iCell] * referenceSpeciesEquation;
    }

    delete []formationRates;
    delete []forwardRates;
    delete []backwardRates;
}

void MixedGas::EnergySource(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble **energySourceTerms, int nCellNumber)
{
    using namespace IDX;
    int nEquation = nl + ntmodel;
    RDouble temperature[3];
    //! The following variable is used to store the source items of species transport equations.
    RDouble *omega = new RDouble[numberOfSpecies];
    RDouble *primitiveVars = new RDouble[nEquation];

    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble referenceEnergyEquation = referenceSpeciesEquation / pow(referenceVelocity, 2);

    //! Compute source term of energy equations.
    RDouble vibrationSource = 0.0, electronSource = 0.0;
    for (int iCell = 0; iCell < nCellNumber; ++ iCell)
    {
        //! Obtain the temperature.
        for (int it = 0; it < ntmodel; ++ it)
        {
            temperature[it] = temperatures[it][iCell];
        }
        if (ntmodel == 1)
        {
            temperature[ITV] = temperature[ITT];
            temperature[ITE] = temperature[ITT];
        }
        else if (ntmodel == 2)
        {
            temperature[ITE] = temperature[ITV];
        }

        //! Obtain the primary variables.
        for (int s = 0; s < nEquation; ++ s)
        {
            primitiveVars[s] = primitiveVariables[s][iCell];
        }
        //! Translate the chemical source term to the dimensional value.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            omega[s] = chemicalSourceTerms[s][iCell] / (cellVolume[iCell] * referenceSpeciesEquation);
        }

        //! Compute source term of the vibration energy equation, and the variables in prim are non-dimensional value.
        if (ntmodel == 2)
        {
            vibrationSource = ComputeVibrationAndElectronEnergySourceTerm(primitiveVars, omega, temperature);
            energySourceTerms[0][iCell] = referenceEnergyEquation * vibrationSource * cellVolume[iCell];
        }
        else if (ntmodel == 3)
        {
            vibrationSource = ComputeVibrationEnergySourceTerm(primitiveVars, omega, temperature);
            electronSource  = ComputeElectronEnergySourceTerm(primitiveVars, omega, temperature);
            energySourceTerms[0][iCell] = referenceEnergyEquation * vibrationSource * cellVolume[iCell];
            energySourceTerms[1][iCell] = referenceEnergyEquation * electronSource * cellVolume[iCell];
        }
    }

    delete []omega;
    delete []primitiveVars;
}

RDouble MixedGas::ComputeVibrationAndElectronEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature)
{
    using namespace IDX;
    //! Obtain the reference values.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble referenceSquareVelocity = referenceVelocity * referenceVelocity;
    //! Obtain the dimensional value of temperature.
    RDouble Ttr = temperature[ITT] * referenceTemperature;
    RDouble Tve = temperature[ITV] * referenceTemperature;
    RDouble density = primitiveVariables[IR] * referenceDensity;
    RDouble pressure = primitiveVariables[IP] * referenceDensity * referenceSquareVelocity;

    RDouble *massFraction = new RDouble[numberOfSpecies];
    RDouble *electronEnergy = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergy = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergyStar1 = new RDouble[numberOfSpecies];
    RDouble *transVibrationTime = new RDouble[numberOfSpecies];
    RDouble *collisionArea = new RDouble[numberOfSpecies];
    //! To obtain the mass fraction of species.
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        massFraction[i] = primitiveVariables[nm + i];
    }
    //! Compute the vibration energy.
    this->ComputeDimensionalVibrationEnergy(temperature[ITV], vibrationEnergy);
    //! Compute the vibration energy computed by translation-rotation temperature.
    this->ComputeDimensionalVibrationEnergy(temperature[ITT], vibrationEnergyStar1);
    //! Compute the electron energy.
    this->ComputeDimensionalElectronEnergy(temperature[ITV], electronEnergy);

    //! nTypes denotes the types of species which is used to identiy the electron, monatom and polyatom.
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();
    //! Obtain the molecular weights of species.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! Compute the relaxation time of translation and vibration.
    GetTranslationAndVibrationRelaxationTime(thermodynamicProperties, massFraction, speciesMass, Ttr, pressure, density, transVibrationTime);

    //! Compute the electron pressure.
    RDouble gasConstant = rjmk;
    RDouble electronPressure = 0.0, electronDensity = 0.0, Ne = 0.0;
    if (this->nElectronIndex >= 0)
    {
        electronDensity = density * massFraction[nElectronIndex];     //! Density of electron.
        //! Compute the number density of electron.
        Ne = AVOGADRO_CONSTANT * electronDensity / speciesMass[nElectronIndex];
        //! Compute the pressure of electron.
        electronPressure = electronDensity * Tve * gasConstant / speciesMass[nElectronIndex];
    }
    //! Compute the effective collision cross-sectional areas of each heavy species.
    GetEffectiveCollisionCrossSectionalArea(thermodynamicProperties, Tve, Ne, collisionArea);

    //! Compute source term of energy term.
    RDouble vibrationOmega = 0.0, fes = 0.0, Ns = 0.0, coef = 0.0, speciesDensity = 0.0;
    RDouble heatVelocity = sqrt(8 * gasConstant * Tve / (PI * this->electronMassDimension));
    coef = 3.0 * gasConstant * electronDensity * (Ttr - Tve);

    for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)   //! Ignore electron.
    {
        speciesDensity = density * massFraction[s];
        Ns = AVOGADRO_CONSTANT * speciesDensity / speciesMass[s];    //! Number density of species.
        fes = Ns * collisionArea[s] * heatVelocity;       //! Frequency of collision.

        //! Compute exchange term of translation-electron energy.
        vibrationOmega += coef * fes / speciesMass[s];

        if (thermodynamicProperties[s].nType < 2)      //! Ignore mono-atom species.
        {
            continue;
        }

        //! Compute exchange term of translation-vibration energy.
        vibrationOmega += speciesDensity * (vibrationEnergyStar1[s] - vibrationEnergy[s]) / transVibrationTime[s];

        //! Compute couple term of chemical reaction and vibration-electron energy.
        vibrationOmega += omega[s] * (vibrationEnergy[s] + electronEnergy[s]);
    }

    delete []massFraction;
    delete []electronEnergy;
    delete []vibrationEnergy;
    delete []vibrationEnergyStar1;
    delete []transVibrationTime;
    delete []collisionArea;
    return vibrationOmega;       //! Return the source item of vibration energy equation.
}

RDouble MixedGas::ComputeVibrationEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature)
{
    using namespace IDX;
    //! Obtain the reference values.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble referenceSquareVelocity = referenceVelocity * referenceVelocity;
    RDouble Ttr = temperature[ITT] * referenceTemperature;
    RDouble Tv =  temperature[ITV] * referenceTemperature;
    RDouble Te =  temperature[ITE] * referenceTemperature;
    RDouble density = primitiveVariables[IR] * referenceDensity;
    RDouble pressure = primitiveVariables[IP] * referenceDensity * referenceSquareVelocity;

    RDouble *vibrationEnergy = new RDouble[numberOfSpecies];
    RDouble *energyStar1 = new RDouble[numberOfSpecies];
    RDouble *energyStar2 = new RDouble[numberOfSpecies];

    RDouble *massFraction = new RDouble[numberOfSpecies];
    RDouble *tvs = new RDouble[numberOfSpecies];
    RDouble *tes = new RDouble[numberOfSpecies];
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        massFraction[s] = primitiveVariables[nm + s];
    }

    //! Compute the vibration energy. This is non-dimensional value.
    this->ComputeDimensionalVibrationEnergy(temperature[ITV], vibrationEnergy);
    //! The vibration energies computed by translation-rotation temperature and electron temperature, respectively.
    this->ComputeDimensionalVibrationEnergy(temperature[ITT], energyStar1);
    this->ComputeDimensionalVibrationEnergy(temperature[ITE], energyStar2);

    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();
    //! Obtain the molecular weights of species.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();

    //! Compute the relaxation time of translation and vibration.
    GetTranslationAndVibrationRelaxationTime(thermodynamicProperties, massFraction, speciesMass, Ttr, pressure, density, tvs);

    //! Compute the electron pressure.
    RDouble gasConstant = rjmk;
    RDouble electronPressure = 0.0, electronDensity = 0.0;
    if (this->nElectronIndex >= 0)
    {
        electronDensity = density * massFraction[nElectronIndex];
        electronPressure = electronDensity * gasConstant * Te / speciesMass[nElectronIndex];
    }
    //! Compute the relaxation time of electron and vibration.
    GetElectronAndVibrationRelaxationTime(thermodynamicProperties, electronPressure, Te, tes);

    RDouble vibrationOmega = 0.0;
    for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
    {
        if (thermodynamicProperties[s].nType < 2)
        {
            continue;
        }
        //! The transport item of translation and vibration energy.
        vibrationOmega += density * massFraction[s] * (energyStar1[s] - vibrationEnergy[s]) / tvs[s];

        //! The coupling term of chemical reaction and vibration energy.
        vibrationOmega += omega[s] * vibrationEnergy[s];
    }
    //! The transport item of vibration and electron energy.Only consider the electric-vibrational effects of Nitrogen.
    if (nNitrogenIndex >= 0)
    {
        vibrationOmega += density * massFraction[nNitrogenIndex] * (energyStar2[nNitrogenIndex] - vibrationEnergy[nNitrogenIndex]) / tes[nNitrogenIndex];
    }

    delete []tvs;
    delete []tes;
    delete []energyStar1;
    delete []energyStar2;
    delete []vibrationEnergy;
    delete []massFraction;
    return vibrationOmega;     //! Return the source item of vibration energy equation.
}

void MixedGas::GetTranslationAndVibrationRelaxationTime(Thermo_Param *thermodynamicProperties, RDouble *massFraction,
                                                        RDouble *speciesMass, RDouble temperature, RDouble pressure, RDouble density, RDouble *tvs)
{
    RDouble tvsMW, tvsP;
    RDouble NA = AVOGADRO_CONSTANT;      //! Avogadro constant(mol-1).
    RDouble P0 = STANDARD_PRESSURE;      //! The standard pressure (Pa).

    RDouble nSum = 0.0;
    RDouble *particle = new RDouble[numberOfSpecies];
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        //! Compute number density of particle.
        particle[i] = density * massFraction[i] * NA / speciesMass[i];
        if (thermodynamicProperties[i].nType > 0)
        {
            nSum += particle[i];
        }
    }

    RDouble coefAsi, coefBsi, coefMsi;
    for (int s = 0; s < this->numberOfSpecies; ++ s)
    {
        if (thermodynamicProperties[s].nType < 1)
        {
            tvs[s] = 1.0e-30;
            continue;
        }

        RDouble fSum = 0.0;
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType < 1)
            {
                continue;
            }
            coefMsi = 1000.0 * (speciesMass[i] * speciesMass[s]) / (speciesMass[i] + speciesMass[s]);      //
            coefAsi = 1.16e-3 * sqrt(coefMsi) * pow(thermodynamicProperties[s].Tvs.Tv[0], 4.0 / 3.0);
            coefBsi = 1.50e-2 * sqrt(sqrt(static_cast<float>(coefMsi)));
            fSum += particle[i] * exp(coefAsi * (pow(temperature, -1.0 / 3.0) - coefBsi) - 18.42);
        }
        tvsMW = (P0 / pressure) * (fSum / nSum);
        tvsP = particle[s] * 3.0e-21 * pow(50000.0 / temperature, 2) * sqrt(8.0 * rjmk * temperature / (PI * speciesMass[s]));
        tvsP = 1.0 / tvsP;

        tvs[s] = tvsMW + tvsP;
    }

    delete []particle;
}

void MixedGas::GetTranslationAndVibrationRelaxationTimeR(Thermo_Param *thermodynamicProperties, RDouble *massFraction,
                                                        RDouble *speciesMass, RDouble temperature, RDouble pressure, RDouble density, RDouble *tvs)
{
    RDouble tvsMW, tvsP;
    RDouble NA = AVOGADRO_CONSTANT;      //! Avogadro constant(mol-1).
    RDouble P0 = STANDARD_PRESSURE;      //! The standard pressure (Pa).

    RDouble nSum = 0.0;
    RDouble *particle = speciesConductivity;   //be carefull
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        //! Compute number density of particle.
        particle[i] = density * massFraction[i] * NA / speciesMass[i];
        if (thermodynamicProperties[i].nType > 0)
        {
            nSum += particle[i];
        }
    }

    RDouble coefAsi, coefBsi, coefMsi;
    for (int s = 0; s < this->numberOfSpecies; ++ s)
    {
        if (thermodynamicProperties[s].nType < 1)
        {
            tvs[s] = 1.0e-30;
            continue;
        }

        RDouble fSum = 0.0;
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType < 1)
            {
                continue;
            }
            coefMsi = 1000.0 * (speciesMass[i] * speciesMass[s]) / (speciesMass[i] + speciesMass[s]);
            coefAsi = 1.16e-3 * sqrt(static_cast<RealData>(coefMsi)) * pow(static_cast<RealData>(thermodynamicProperties[s].Tvs.Tv[0]), static_cast<RealData>(4.0 / 3.0));
            coefBsi = 1.50e-2 * sqrt(sqrt(static_cast<RealData>(coefMsi)));
            fSum += particle[i] * exp(static_cast<RealData>(coefAsi * (pow(static_cast<RealData>(temperature), static_cast<RealData>(-1.0 / 3.0)) - coefBsi) - 18.42));
        }
        tvsMW = (P0 / pressure) * (fSum / nSum);
        tvsP = 50000.0 / temperature;
        tvsP = particle[s] * 3.0e-21 * tvsP* tvsP * sqrt(static_cast<RealData>(8.0 * rjmk * temperature / (PI * speciesMass[s])));
        tvsP = 1.0 / tvsP;

        tvs[s] = tvsMW + tvsP;
    }
}

void MixedGas::GetElectronAndVibrationRelaxationTime(Thermo_Param *thermodynamicProperties, RDouble electronPressure, RDouble electronTemperature, RDouble *tes)
{
    if (nNitrogenIndex < 0)
    {
        return;
    }
    //! Only consider the electric-vibrational effects of Nitrogen.
    RDouble kB = BOLTZMANN_CONSTANT;      //! Boltzmann constant.
    RDouble ce = ELECTRIC_QUANTITY;       //! The electric quantity.

    //! Compute the vibration excitation coefficient.
    RDouble fx = kB * electronTemperature / ce;
    //! Obtain the vibration excitation coefficient.
    RDouble kevj, kevj10 = 0.0;

    RDouble *av = this->CurveFits->GetNitrogenVibrationalExcitationRateA();
    RDouble *bv = this->CurveFits->GetNitrogenVibrationalExcitationRateB();
    RDouble *cv = this->CurveFits->GetNitrogenVibrationalExcitationRateC();

    for (int i = 0; i < 10; ++ i)
    {
        kevj = 1.0e-15 * av[i] * pow(fx, -1.5) * exp(bv[i] / fx + cv[i]);
        kevj10 += kevj * pow(static_cast <RDouble> (i + 1), 2);
    }

    //! Compute the electric-vibrational effects of Nitrogen.
    tes[nNitrogenIndex] = kB * electronTemperature / (electronPressure + SMALL);
    tes[nNitrogenIndex] /= 0.5 * kevj10 * (1.0 - exp(-thermodynamicProperties[nNitrogenIndex].Tvs.Tv[0] / electronTemperature));
}

RDouble MixedGas::ComputeElectronEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature)
{
    using namespace IDX;
    //! Obtain the reference values.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble referenceSquareVelocity = referenceVelocity * referenceVelocity;
    RDouble Ttr = temperature[ITT] * referenceTemperature;
    RDouble Tv  = temperature[ITV] * referenceTemperature;
    RDouble Te  = temperature[ITE] * referenceTemperature;
    RDouble density = primitiveVariables[IR] * referenceDensity;
    RDouble pressure = primitiveVariables[IP] * referenceDensity * referenceSquareVelocity;

    RDouble *energyStar2 = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergy = new RDouble[numberOfSpecies];
    RDouble *electronEnergy = new RDouble[numberOfSpecies];
    RDouble *massFraction = new RDouble[numberOfSpecies];
    RDouble *tes = new RDouble[numberOfSpecies];
    RDouble *sigma = new RDouble[numberOfSpecies];
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        massFraction[s] = primitiveVariables[nm + s];
    }

    //! Obtain the vibration energy computed by the electron temperature.
    this->ComputeDimensionalVibrationEnergy(temperature[ITE], energyStar2);
    //! Compute the vibration and electron energy. Both are the non-dimensional value.
    this->ComputeDimensionalVibrationEnergy(temperature[ITV], vibrationEnergy);
    //! Compute the electron energy.
    this->ComputeDimensionalElectronEnergy(temperature[ITE], electronEnergy);

    //! nTypes denotes the types of species which is used to identiy the electron, monatom and polyatom.
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();

    // Compute the electron pressure.
    RDouble gasConstant = rjmk;
    RDouble electronPressure = 0.0, Ne = 0.0, electronDensity = 0.0;
    if (this->nElectronIndex >= 0)
    {
        electronDensity = density * massFraction[nElectronIndex];
        //! Compute the number density of electron.
        Ne = AVOGADRO_CONSTANT * electronDensity / speciesMass[nElectronIndex];
        //! Compute the pressure of electron.
        electronPressure = electronDensity * Te * gasConstant / speciesMass[nElectronIndex];
    }

    //! Compute the relaxation time of electron and vibration.
    GetElectronAndVibrationRelaxationTime(thermodynamicProperties, electronPressure, Te, tes);
    //! Compute the effective collision cross-sectional areas of each heavy species.
    GetEffectiveCollisionCrossSectionalArea(thermodynamicProperties, Te, Ne, sigma);

    RDouble electronOmega = 0.0, fes = 0.0, Ns = 0.0, speciesDensity = 0.0;
    RDouble heatVelocity = sqrt(8.0 * gasConstant * Te / (PI * this->electronMassDimension));
    RDouble coef = 3.0 * electronDensity * gasConstant * (Ttr - Te);

    for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
    {
        speciesDensity = density * massFraction[s];
        //! Compute the number density of the current species.
        Ns = AVOGADRO_CONSTANT * speciesDensity / speciesMass[s];
        //! Compute the frequence of effective collision between the electron and s-th heavy particle.
        fes = Ns * sigma[s] * heatVelocity; //frequency of collision.

        //! Compute the exchange term of translation and electron energy.
        electronOmega += coef * fes / speciesMass[s];

        //! Compute the couple term of chemical reaction and electron energy.
        electronOmega += omega[s] * electronEnergy[s];
    }
    //! The transport item of vibration and electron energy.Only consider the electric-vibrational effects of Nitrogen.
    if (nNitrogenIndex >= 0)
    {
        electronOmega -= density * massFraction[nNitrogenIndex] * (energyStar2[nNitrogenIndex] - vibrationEnergy[nNitrogenIndex]) / tes[nNitrogenIndex];
    }
    delete []energyStar2;
    delete []vibrationEnergy;
    delete []electronEnergy;
    delete []massFraction;
    delete []tes;
    delete []sigma;
    return electronOmega;     //! Return the source item of electron energy equation.
}

void MixedGas::GetEffectiveCollisionCrossSectionalArea(Thermo_Param *thermodynamicProperties, RDouble electronTemperature, RDouble electronNumberDensity, RDouble *sigma)
{
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    RDouble ce = ELECTRIC_QUANTITY;       //! Obtain the electric quantity.
    RDouble kB = BOLTZMANN_CONSTANT;      //! Obtain the Boltzmann constant. 1 J/K = 1 C*V/K.
    RDouble eps = VACUUM_PERMITTIVITY;    //! Obtain the permittivity of vacuum.

    RDouble *aes = CurveFits->GetEffectiveCollisionCrossSectionalAreaCurveFitsA();
    RDouble *bes = CurveFits->GetEffectiveCollisionCrossSectionalAreaCurveFitsB();
    RDouble *ces = CurveFits->GetEffectiveCollisionCrossSectionalAreaCurveFitsC();
    RDouble T1 = electronTemperature;
    RDouble T2 = T1 * T1;

    for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)    //! Ignore the electron.
    {
        if (ionType[s] == 0)        //! The neutral particle.
        {
            sigma[s] = aes[s] + bes[s] * T1 + ces[s] * T2;
        }

        if (ionType[s] > 0)      //! The positive particle or ion.
        {
            RDouble tmp = (ce * ce) / (eps * kB * electronTemperature);
            //! Compute the scattering factor of 90-degree angle.
            RDouble b0 = ionType[s] * tmp / (12.0 * PI);
            //! Compute the radius of Debye.
            RDouble rD = sqrt(1.0 / tmp / (electronNumberDensity + SMALL));

            sigma[s] = (8.0 * PI / 3.0) * (b0 * b0) * log(1.0 + pow(rD / b0, 2));
        }
    }
}

//void RealGas::ChemicalSpectrumRadius(RDouble **q, RDouble **t, RDouble *vol, RDouble **srs, int nlen)
//{
//    using namespace IDX;
//    RDouble *molecularWeightDimensional  = molecularProperty->GetMolecularWeightDimensional();
//    RDouble *oMolecularWeightDimensional = molecularProperty->GetMolecularWeightReciprocalDimensional();
//
//    RDouble refDimensionalTemperature = this->referenceParameterDimensional.GetTemperature();
//    RDouble reference_velocity_dimensional = this->referenceParameterDimensional.GetVelocity();
//    RDouble refDimensionalDensity = this->referenceParameterDimensional.GetDensity();
//
//    RDouble temperatureDimensional, densityDimensional, oTemperatureDimensional, oTemperatureDimensional2;
//    RDouble diffp, diffn, posr, negr, comm, pztt;
//    RDouble xis, srst;
//    RDouble density, temperature;
//    RDouble odensityDimensionalJS, odensityDimensionalNS;
//    RDouble cpDimensional, cvDimensional, omav, coef_tr;
//
//    RDouble lov = reynoldsReferenceLengthDimensional / reference_velocity_dimensional;
//
//    int neqn = nl + nchem;
//
//    RDouble *afr = reactionRate->GetForwardReactionRateCoefficientA();
//    RDouble *bfr = reactionRate->GetForwardReactionRateCoefficientB();
//    RDouble *cfr = reactionRate->GetForwardReactionRateCoefficientC();
//             
//    RDouble *abr = reactionRate->GetBackwardReactionRateCoefficientA();
//    RDouble *bbr = reactionRate->GetBackwardReactionRateCoefficientB();
//    RDouble *cbr = reactionRate->GetBackwardReactionRateCoefficientC();
//
//    int **cvp = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
//    int **cvn = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
//    RDouble **cpz = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();
//
//    RDouble *cq = new RDouble [numberOfReaction];
//    RDouble *kf = new RDouble [numberOfReaction];
//    RDouble *kb = new RDouble [numberOfReaction];
//
//    RDouble *pos = new RDouble[numberOfReaction];
//    RDouble *neg = new RDouble[numberOfReaction];
//    RDouble *pzt = new RDouble[numberOfReaction];
//    int    *npz = new int   [numberOfReaction];
//
//    RDouble *c_a1 = new RDouble[numberOfSpecies];
//    RDouble *c_a2 = new RDouble[numberOfSpecies];
//    RDouble *c_bk = new RDouble[numberOfSpecies];
//    RDouble *c_tr = new RDouble[numberOfSpecies];
//
//    RDouble *prim = new RDouble[neqn];
//
//    for (int i = 0; i < nlen; ++ i)
//    {
//        density     = q[IR ][i];
//        temperature = t[ITT][i];
//
//        densityDimensional = density * refDimensionalDensity;
//        temperatureDimensional = temperature * refDimensionalTemperature;
//        oTemperatureDimensional = one / ( temperatureDimensional + SMALL );
//        oTemperatureDimensional2 = SQR( oTemperatureDimensional );
//
//        for (int m = 0; m < neqn; ++ m)
//        {
//            prim[m] = q[m][i];
//        }
//
//        for (int ireaction = 0; ireaction < numberOfReaction; ++ ireaction)
//        {
//            pos[ireaction] = 1.0;
//            neg[ireaction] = 1.0;
//            pzt[ireaction] = 0.0;
//            for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
//            {
//                xis = densityDimensional * ( prim[nm+ispecies] + SMALL ) * oMolecularWeightDimensional[ispecies];
//
//                for (int nexp = 0; nexp < cvp[ireaction][ispecies]; ++ nexp)
//                {
//                    pos[ireaction] *= xis;
//                }
//
//                for (int nexp = 0; nexp < cvn[ireaction][ispecies]; ++ nexp)
//                {
//                    neg[ireaction] *= xis;
//                }
//
//                pzt[ireaction] += xis * cpz[ireaction][ispecies];
//            }
//
//            npz[ireaction] = 0;
//            for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
//            {
//                npz[ireaction] += static_cast<int>(cpz[ireaction][ispecies] + 0.1);
//            }
//
//            if (npz[ireaction] == 0)
//            {
//                pzt[ireaction] = one;
//            }
//
//            kf[ireaction] = afr[ireaction] * pow(temperatureDimensional, bfr[ireaction]) * exp(-cfr[ireaction] * oTemperatureDimensional);
//            kb[ireaction] = abr[ireaction] * pow(temperatureDimensional, bbr[ireaction]) * exp(-cbr[ireaction] * oTemperatureDimensional);
//        }
//
//        for (int ispecies = 0; ispecies < numberOfSpecies - 1; ++ ispecies)
//        {
//            c_a2[ispecies] = 0.0;
//            for (int ireaction = 0; ireaction < numberOfReaction; ++ ireaction)
//            {
//                diffp = bfr[ireaction] * oTemperatureDimensional + cfr[ireaction] * oTemperatureDimensional2;
//                diffn = bbr[ireaction] * oTemperatureDimensional + cbr[ireaction] * oTemperatureDimensional2;
//                posr = pos[ireaction] * pzt[ireaction] * kf[ireaction] * diffp;
//                negr = neg[ireaction] * pzt[ireaction] * kb[ireaction] * diffn;
//                comm = ( cvn[ireaction][ispecies] - cvp[ireaction][ispecies] );
//                c_a2[ispecies] += comm * ( posr - negr );
//            }
//            c_a2[ispecies] *= molecularWeightDimensional[ispecies];
//        }
//
//        ComputeSpeciesEnthalpyDimensional(temperature, speciesEnthalpy);
//        ComputeSpeciesConstantPressureSpecificHeatDimensional(temperature, speciesCp);
//        ComputeMixtureByPrimitive(prim, speciesCp, cpDimensional);
//        ComputeMixtureByPrimitive(prim, oMolecularWeightDimensional, omav);
//
//        cvDimensional = cpDimensional - rjmk * omav;
//        coef_tr       = - 1.0 / ( densityDimensional * cvDimensional );
//
//        for (int ispecies = 0; ispecies < numberOfSpecies - 1; ++ ispecies)
//        {
//            c_tr[ispecies] = coef_tr * ( speciesEnthalpy[ispecies] - speciesEnthalpy[numberOfSpecies-1] - rjmk * temperatureDimensional * 
//                                 ( oMolecularWeightDimensional[ispecies] - oMolecularWeightDimensional[numberOfSpecies-1] )  );
//        }
//
//        for (int ispecies = 0; ispecies < numberOfSpecies - 1; ++ ispecies)
//        {
//            for (int jspecies = 0; jspecies < numberOfSpecies - 1; ++ jspecies)
//            {
//                c_bk[jspecies] = 0.0;
//                odensityDimensionalJS = 1.0 / ( densityDimensional * ( prim[nm+jspecies         ] + SMALL ) );
//                odensityDimensionalNS = 1.0 / ( densityDimensional * ( prim[nm+numberOfSpecies-1] + SMALL ) );
//                for (int ireaction = 0; ireaction < numberOfReaction; ++ ireaction)
//                {
//                    diffp = cvp[ireaction][jspecies] * odensityDimensionalJS - cvp[ireaction][numberOfSpecies-1] * odensityDimensionalNS;
//                    diffn = cvn[ireaction][jspecies] * odensityDimensionalJS - cvn[ireaction][numberOfSpecies-1] * odensityDimensionalNS;
//                    posr  = pos[ireaction] * kf[ireaction] * diffp;
//                    negr  = neg[ireaction] * kb[ireaction] * diffn;
//                    comm  = ( cvn[ireaction][ispecies] - cvp[ireaction][ispecies] );
//                    c_bk[jspecies] += comm * ( posr - negr ) * pzt[ireaction];
//                    if (npz[ireaction] != 0)
//                    {
//                        pztt = cpz[ireaction][jspecies         ] * oMolecularWeightDimensional[jspecies         ]
//                             - cpz[ireaction][numberOfSpecies-1] * oMolecularWeightDimensional[numberOfSpecies-1];
//                        c_bk[jspecies] += comm * ( kf[ireaction] * pos[ireaction] - kb[ireaction] * neg[ireaction] ) * pztt;
//                    }
//                }
//                c_bk[jspecies] = molecularWeightDimensional[ispecies] * c_bk[jspecies];
//            }
//            srst = 0.0;
//            for (int jspecies = 0; jspecies < numberOfSpecies - 1; ++ jspecies)
//            {
//                srst += ABS( c_bk[jspecies] + c_a2[ispecies] * c_tr[jspecies] );
//            }
//            srs[ispecies][i] = srst * lov * vol[i];
//        }
//    }
//    delete []cq;
//    delete []kf;
//    delete []kb;
//
//    delete []pos;
//    delete []neg;
//    delete []pzt;
//    delete []npz;
//
//    delete []c_a1;
//    delete []c_a2;
//    delete []c_bk;
//    delete []c_tr;
//
//    delete []prim;
//}

void MixedGas::ChemicalSpectrumRadius(RDouble **q, RDouble **t, RDouble *vol, RDouble **srs, int iCell)
{
    using namespace IDX;
    RDouble *molecularWeightDimensional  = molecularProperty->GetMolecularWeightDimensional();
    RDouble *oMolecularWeightDimensional = molecularProperty->GetMolecularWeightReciprocalDimensional();

    RDouble refDimensionalTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble reference_velocity_dimensional = this->referenceParameterDimensional.GetVelocity();
    RDouble refDimensionalDensity = this->referenceParameterDimensional.GetDensity();

    RDouble temperatureDimensional, densityDimensional, oTemperatureDimensional, oTemperatureDimensional2;
    RDouble diffp, diffn, posr, negr, comm, pztt;
    RDouble xis, srst;
    RDouble density, temperature;
    RDouble odensityDimensionalJS, odensityDimensionalNS;
    RDouble cpDimensional, cvDimensional, omav, coef_tr;

    RDouble lov = reynoldsReferenceLengthDimensional / reference_velocity_dimensional;

    int neqn = nl + nchem;

    RDouble *afr = reactionRate->GetForwardReactionRateCoefficientA();
    RDouble *bfr = reactionRate->GetForwardReactionRateCoefficientB();
    RDouble *cfr = reactionRate->GetForwardReactionRateCoefficientC();

    RDouble *abr = reactionRate->GetBackwardReactionRateCoefficientA();
    RDouble *bbr = reactionRate->GetBackwardReactionRateCoefficientB();
    RDouble *cbr = reactionRate->GetBackwardReactionRateCoefficientC();

    int **cvp = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
    int **cvn = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
    RDouble **cpz = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();

    RDouble *cq = new RDouble [numberOfReaction];
    RDouble *kf = new RDouble [numberOfReaction];
    RDouble *kb = new RDouble [numberOfReaction];

    RDouble *pos = new RDouble[numberOfReaction];
    RDouble *neg = new RDouble[numberOfReaction];
    RDouble *pzt = new RDouble[numberOfReaction];
    int    *npz = new int   [numberOfReaction];

    RDouble *c_a1 = new RDouble[numberOfSpecies];
    RDouble *c_a2 = new RDouble[numberOfSpecies];
    RDouble *c_bk = new RDouble[numberOfSpecies];
    RDouble *c_tr = new RDouble[numberOfSpecies];

    RDouble *prim = new RDouble[neqn];

    density     = q[IR ][iCell];
    temperature = t[ITT][iCell];

    densityDimensional = density * refDimensionalDensity;
    temperatureDimensional = temperature * refDimensionalTemperature;
    oTemperatureDimensional = one / ( temperatureDimensional + SMALL );
    oTemperatureDimensional2 = SQR( oTemperatureDimensional );

    for (int m = 0; m < neqn; ++ m)
    {
        prim[m] = q[m][iCell];
    }

    for (int ireaction = 0; ireaction < numberOfReaction; ++ ireaction)
    {
        pos[ireaction] = 1.0;
        neg[ireaction] = 1.0;
        pzt[ireaction] = 0.0;
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            xis = densityDimensional * ( prim[nm+ispecies] + SMALL ) * oMolecularWeightDimensional[ispecies];

            for (int nexp = 0; nexp < cvp[ireaction][ispecies]; ++ nexp)
            {
                pos[ireaction] *= xis;
            }

            for (int nexp = 0; nexp < cvn[ireaction][ispecies]; ++ nexp)
            {
                neg[ireaction] *= xis;
            }

            pzt[ireaction] += xis * cpz[ireaction][ispecies];
        }

        npz[ireaction] = 0;
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            npz[ireaction] += static_cast<int>(cpz[ireaction][ispecies] + 0.1);
        }

        if (npz[ireaction] == 0)
        {
            pzt[ireaction] = one;
        }

        kf[ireaction] = afr[ireaction] * pow(temperatureDimensional, bfr[ireaction]) * exp(-cfr[ireaction] * oTemperatureDimensional);
        kb[ireaction] = abr[ireaction] * pow(temperatureDimensional, bbr[ireaction]) * exp(-cbr[ireaction] * oTemperatureDimensional);
    }

    for (int ispecies = 0; ispecies < numberOfSpecies - 1; ++ ispecies)
    {
        c_a2[ispecies] = 0.0;
        for (int ireaction = 0; ireaction < numberOfReaction; ++ ireaction)
        {
            diffp = bfr[ireaction] * oTemperatureDimensional + cfr[ireaction] * oTemperatureDimensional2;
            diffn = bbr[ireaction] * oTemperatureDimensional + cbr[ireaction] * oTemperatureDimensional2;
            posr = pos[ireaction] * pzt[ireaction] * kf[ireaction] * diffp;
            negr = neg[ireaction] * pzt[ireaction] * kb[ireaction] * diffn;
            comm = ( cvn[ireaction][ispecies] - cvp[ireaction][ispecies] );
            c_a2[ispecies] += comm * ( posr - negr );
        }
        c_a2[ispecies] *= molecularWeightDimensional[ispecies];
    }

    ComputeSpeciesEnthalpyDimensional(temperature, speciesEnthalpy);
    ComputeSpeciesConstantPressureSpecificHeatDimensional(temperature, speciesCp);
    ComputeMixtureByPrimitive(prim, speciesCp, cpDimensional);
    ComputeMixtureByPrimitive(prim, oMolecularWeightDimensional, omav);

    cvDimensional = cpDimensional - rjmk * omav;
    coef_tr       = - 1.0 / ( densityDimensional * cvDimensional );

    for (int ispecies = 0; ispecies < numberOfSpecies - 1; ++ ispecies)
    {
        c_tr[ispecies] = coef_tr * ( speciesEnthalpy[ispecies] - speciesEnthalpy[numberOfSpecies-1] - rjmk * temperatureDimensional * 
            ( oMolecularWeightDimensional[ispecies] - oMolecularWeightDimensional[numberOfSpecies-1] )  );
    }

    for (int ispecies = 0; ispecies < numberOfSpecies - 1; ++ ispecies)
    {
        for (int jspecies = 0; jspecies < numberOfSpecies - 1; ++ jspecies)
        {
            c_bk[jspecies] = 0.0;
            odensityDimensionalJS = 1.0 / ( densityDimensional * ( prim[nm+jspecies         ] + SMALL ) );
            odensityDimensionalNS = 1.0 / ( densityDimensional * ( prim[nm+numberOfSpecies-1] + SMALL ) );
            for (int ireaction = 0; ireaction < numberOfReaction; ++ ireaction)
            {
                diffp = cvp[ireaction][jspecies] * odensityDimensionalJS - cvp[ireaction][numberOfSpecies-1] * odensityDimensionalNS;
                diffn = cvn[ireaction][jspecies] * odensityDimensionalJS - cvn[ireaction][numberOfSpecies-1] * odensityDimensionalNS;
                posr  = pos[ireaction] * kf[ireaction] * diffp;
                negr  = neg[ireaction] * kb[ireaction] * diffn;
                comm  = ( cvn[ireaction][ispecies] - cvp[ireaction][ispecies] );
                c_bk[jspecies] += comm * ( posr - negr ) * pzt[ireaction];
                if (npz[ireaction] != 0)
                {
                    pztt = cpz[ireaction][jspecies         ] * oMolecularWeightDimensional[jspecies         ]
                    - cpz[ireaction][numberOfSpecies-1] * oMolecularWeightDimensional[numberOfSpecies-1];
                    c_bk[jspecies] += comm * ( kf[ireaction] * pos[ireaction] - kb[ireaction] * neg[ireaction] ) * pztt;
                }
            }
            c_bk[jspecies] = molecularWeightDimensional[ispecies] * c_bk[jspecies];
        }
        srst = 0.0;
        for (int jspecies = 0; jspecies < numberOfSpecies - 1; ++ jspecies)
        {
            srst += ABS( c_bk[jspecies] + c_a2[ispecies] * c_tr[jspecies] );
        }
        srs[ispecies][iCell] = srst * lov * vol[iCell];
    }

    delete []cq;
    delete []kf;
    delete []kb;

    delete []pos;
    delete []neg;
    delete []pzt;
    delete []npz;

    delete []c_a1;
    delete []c_a2;
    delete []c_bk;
    delete []c_tr;

    delete []prim;
}


void MixedGas::ComputeChemicalSourceDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{
    using namespace IDX;     //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();

    //! For source term in species transport equation.
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble reference1 = referenceDensity * referenceSpeciesEquation;
    RDouble reference2 = referenceVelocity * reynoldsReferenceLengthDimensional;

    RDouble *forwardRateCoefA = reactionRate->GetForwardReactionRateCoefficientA();
    RDouble *forwardRateCoefB = reactionRate->GetForwardReactionRateCoefficientB();
    RDouble *forwardRateCoefC = reactionRate->GetForwardReactionRateCoefficientC();
    RDouble *backwardRateCoefA = reactionRate->GetBackwardReactionRateCoefficientA();
    RDouble *backwardRateCoefB = reactionRate->GetBackwardReactionRateCoefficientB();
    RDouble *backwardRateCoefC = reactionRate->GetBackwardReactionRateCoefficientC();

    int **forwardStoichNumber = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
    int **backwardStoichNumber = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
    RDouble **thirdBodyStoichNumber = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();

    //! Record the coefficients in the chemical reactions.
    RDouble *forwardRates = new RDouble [numberOfReaction];
    RDouble *backwardRates = new RDouble [numberOfReaction];
    RDouble *collisionRates = new RDouble [numberOfReaction];
    RDouble *energy = new RDouble [numberOfSpecies];
    RDouble *massFraction = new RDouble [numberOfSpecies];
    RDouble *densityMassFraction = new RDouble [numberOfSpecies];
    int nEquation = this->nm + this->numberOfSpecies + this->ntmodel - 1;
    RDouble **partialValues = NewPointer2 <RDouble> (numberOfReaction, nEquation);

    //! The dimensional value of control temperature of chemical reaction.
    RDouble forwardReactionTemperature, backwardReactionTemperature;
    RDouble densityPerMass, density, mixGasCv, kf, kb, forwardEnergy, backwardEnergy;
    RDouble gamaj = 0.0, Drf = 0.0, Drb = 0.0, phir = 0.0, theta[3] = {0.0}, coefA = 0.0;
    int electronIndex = numberOfSpecies - 1; // The index of electron in the array of species.
    if (nElectronIndex >= 0)
    {
        electronIndex = nElectronIndex;
    }

    //! Compute the elements of Jacobian matrix.
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();
    for (int iCell = 0; iCell < nCellNumber; ++ iCell)
    {
        //! The control temperature of front reaction.
        forwardReactionTemperature = MAX(temperatures[ITT][iCell] * referenceTemperature, SMALL);
        //! The control temperature of back reaction. it is temperature of translation and rotation.
        backwardReactionTemperature = forwardReactionTemperature;
        if (ntmodel > 1)
        {
            forwardReactionTemperature = pow(forwardReactionTemperature, 0.6) * pow(temperatures[ITV][iCell] * referenceTemperature, 0.4);
        }

        density = primitiveVariables[IR][iCell] * referenceDensity;
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            massFraction[s] = primitiveVariables[nm + s][iCell];
            densityMassFraction[s] = MAX(SMALL, density * massFraction[s]);
        }

        if (ntmodel == 1)      //! One-Temperature Model.
        {
            //! Compute the specific heat at constant volume of each species.
            this->ComputeSpeciesConstantVolumeSpecificHeatDimensional(temperatures[ITT][iCell], energy);
            this->ComputeMixtureByMassFraction(massFraction, energy, mixGasCv);
            for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute internal energy of each species.
            {
                energy[s] = energy[s] * forwardReactionTemperature;
            }
        }
        else         //! Multi-Temperature Model.
        {
            // Compute the specific heat at constant volume of each species.
            this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(energy);
            this->ComputeMixtureByMassFraction(massFraction, energy, mixGasCv);
            for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute internal energy of each species.
            {
                energy[s] = energy[s] * forwardReactionTemperature + thermodynamicProperties[s].hs0;
            }
        }
        forwardEnergy = density * mixGasCv * forwardReactionTemperature;
        backwardEnergy = density * mixGasCv * backwardReactionTemperature;

        for (int r = 0; r < numberOfReaction; ++ r)
        {
            forwardRates[r] = 1.0;
            backwardRates[r] = 1.0;
            collisionRates[r] = 0.0;
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                densityPerMass = densityMassFraction[s] * speciesMassReciprocal[s];
                for (int i = 0; i < forwardStoichNumber[r][s]; ++ i)
                {
                    forwardRates[r] *= densityPerMass;
                }
                for (int i = 0; i < backwardStoichNumber[r][s]; ++ i)
                {
                    backwardRates[r] *= densityPerMass;
                }
                collisionRates[r] += densityPerMass * thirdBodyStoichNumber[r][s];
            }
            if (fabs(collisionRates[r]) < SMALL) // RCr
            {
                collisionRates[r] = 1.0;
            }

            //! Compute equilibrium constant of the front reaction.
            kf = forwardRateCoefA[r] * pow(forwardReactionTemperature, forwardRateCoefB[r]) * exp(-forwardRateCoefC[r] / forwardReactionTemperature);
            //! Compute equilibrium constant of the reverse reaction.
            kb = backwardRateCoefA[r] * pow(backwardReactionTemperature, backwardRateCoefB[r]) * exp(-backwardRateCoefC[r] / backwardReactionTemperature);
            forwardRates[r] = kf * forwardRates[r];     //! RAr.
            backwardRates[r] = kb * backwardRates[r];   //! RBr.

            Drf = (forwardRateCoefB[r] + forwardRateCoefC[r] / forwardReactionTemperature) / forwardEnergy;
            Drb = (backwardRateCoefB[r] + backwardRateCoefC[r] / backwardReactionTemperature) / backwardEnergy;

            //! Compute the partial derivatives to each species.
            for (int s = 0; s < this->nLastSpeciesIndex; ++ s)
            {
                //! One-Temperature Model.
                coefA = ionType[s] * speciesMassReciprocal[s] * this->electronMassDimension;
                gamaj = energy[nLastSpeciesIndex] - energy[s] + coefA * (energy[nLastSpeciesIndex] - energy[electronIndex]);

                theta[0]  = forwardStoichNumber[r][s] / densityMassFraction[s] - forwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex];
                theta[0] += coefA * (forwardStoichNumber[r][electronIndex] / densityMassFraction[electronIndex] - forwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex]);
                theta[0] += Drf * gamaj;
                theta[0]  = forwardRates[r] * theta[0];

                theta[1]  = backwardStoichNumber[r][s] / densityMassFraction[s] - backwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex];
                theta[1] += coefA * (backwardStoichNumber[r][electronIndex] / densityMassFraction[electronIndex] - backwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex]);
                theta[1] += Drb * gamaj;
                theta[1]  = backwardRates[r] * theta[1];

                theta[2]  = thirdBodyStoichNumber[r][s] * speciesMassReciprocal[s] - thirdBodyStoichNumber[r][nLastSpeciesIndex] * speciesMassReciprocal[nLastSpeciesIndex];
                theta[2] += coefA * (thirdBodyStoichNumber[r][electronIndex] * speciesMassReciprocal[electronIndex] - thirdBodyStoichNumber[r][nLastSpeciesIndex] * speciesMassReciprocal[nLastSpeciesIndex]);

                //! Compute the partial derivatives to the species j.
                partialValues[r][nm + s] = (theta[0] - theta[1]) * collisionRates[r] + (forwardRates[r] - backwardRates[r]) * theta[2];
            }
            //! Compute the partial derivative to vibration and electron energy.
            for (int it = 1; it < ntmodel; ++ it)
            {
                phir = (Drf * forwardRates[r] - Drb * backwardRates[r]) * collisionRates[r]; // (dRr / dEv).
                partialValues[r][nl + it] = -phir;
            }
        }

        //! Compute the chemical source term.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            chemicalSourceTerms[s][iCell] = 0.0;
            for (int r = 0; r < numberOfReaction; ++ r)
            {
                chemicalSourceTerms[s][iCell] += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * (forwardRates[r] - backwardRates[r]) * collisionRates[r];
            }
            chemicalSourceTerms[s][iCell] *= referenceSpeciesEquation * speciesMass[s] * cellVolume[iCell];
        }

        //! Construct elements of the Jacobian matrix of the chemical source term.
        for (int s = 0; s < this->nLastSpeciesIndex; ++ s)
        {
            sourceDerivatives[s][nm + s][iCell] = 0.0;
            for (int r = 0; r < numberOfReaction; ++ r)
            {
                sourceDerivatives[s][nm + s][iCell] += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * partialValues[r][nm + s]; // (dRr / dcs).
            }
            sourceDerivatives[s][nm + s][iCell] *= reference1 * speciesMass[s] * cellVolume[iCell];// Translate to the nondimensional values.
        }
        for (int it = 1; it < ntmodel; ++ it)      //! Compute the partial derivatives of species source term to the vibration energy.
        {
            for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
            {
                sourceDerivatives[s][nl + it][iCell] = 0.0;
                for (int r = 0; r < numberOfReaction; ++ r)
                {
                    sourceDerivatives[s][nl + it][iCell] += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * partialValues[r][nl + it]; // (dRr / dEv).
                }
                sourceDerivatives[s][nl + it][iCell] *= reference2 * speciesMass[s] * cellVolume[iCell];
            }
        }
    }

    delete []energy;
    delete []forwardRates;
    delete []backwardRates;
    delete []collisionRates;
    delete []massFraction;
    delete []densityMassFraction;
    DelPointer2(partialValues);
}

void MixedGas::ComputeChemicalSourceJacobian(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives)
{
    using namespace IDX;     //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();

    //! For source term in species transport equation.
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble reference1 = referenceDensity * referenceSpeciesEquation;
    RDouble reference2 = referenceVelocity * reynoldsReferenceLengthDimensional;

    RDouble *forwardRateCoefA = reactionRate->GetForwardReactionRateCoefficientA();
    RDouble *forwardRateCoefB = reactionRate->GetForwardReactionRateCoefficientB();
    RDouble *forwardRateCoefC = reactionRate->GetForwardReactionRateCoefficientC();
    RDouble *backwardRateCoefA = reactionRate->GetBackwardReactionRateCoefficientA();
    RDouble *backwardRateCoefB = reactionRate->GetBackwardReactionRateCoefficientB();
    RDouble *backwardRateCoefC = reactionRate->GetBackwardReactionRateCoefficientC();

    int **forwardStoichNumber = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
    int **backwardStoichNumber = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
    RDouble **thirdBodyStoichNumber = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();

    //! Record the coefficients in the chemical reactions.
    RDouble *forwardRates = new RDouble [numberOfReaction];
    RDouble *backwardRates = new RDouble [numberOfReaction];
    RDouble *collisionRates = new RDouble [numberOfReaction];
    RDouble *energy = new RDouble [numberOfSpecies];
    RDouble *densityMassFraction = new RDouble [numberOfSpecies];
    int nEquation = this->nm + this->numberOfSpecies + this->ntmodel - 1;
    RDouble **partialValues = NewPointer2 <RDouble> (numberOfReaction, nEquation);

    //! The dimensional value of control temperature of chemical reaction;
    RDouble forwardReactionTemperature, backwardReactionTemperature;
    RDouble densityPerMass, density, mixGasCv, kf, kb, forwardEnergy, backwardEnergy;
    RDouble gamaj = 0.0, Drf = 0.0, Drb = 0.0, phir = 0.0, theta[3] = {0.0}, coefA = 0.0;
    int electronIndex = numberOfSpecies - 1; // The index of electron in the array of species.
    if (nElectronIndex >= 0)
    {
        electronIndex = nElectronIndex;
    }

    //! Compute the elements of Jacobian matrix.
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();

    //! The control temperature of front reaction.
    forwardReactionTemperature = MAX(temperatures[ITT] * referenceTemperature, SMALL);
    //! The control temperature of back reaction. it is temperature of translation and rotation.
    backwardReactionTemperature = forwardReactionTemperature;
    if (ntmodel > 1)
    {
        forwardReactionTemperature = pow(forwardReactionTemperature, 0.6) * pow(temperatures[ITV] * referenceTemperature, 0.4);
    }

    density = primitiveVariables[IR] * referenceDensity;
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        massFractions[s] = primitiveVariables[nm + s];
        densityMassFraction[s] = MAX(SMALL, density * massFractions[s]);
    }

    if (ntmodel == 1)      //! One-Temperature Model.
    {
        //! Compute the specific heat at constant volume of each species.
        this->ComputeSpeciesConstantVolumeSpecificHeatDimensional(temperatures[ITT], energy);
        this->ComputeMixtureByMassFraction(massFractions, energy, mixGasCv);
        for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute internal energy of each species.
        {
            energy[s] = energy[s] * forwardReactionTemperature;
        }
    }
    else         //! Multi-Temperature Model.
    {
        // Compute the specific heat at constant volume of each species.
        this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(energy);
        this->ComputeMixtureByMassFraction(massFractions, energy, mixGasCv);
        for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute internal energy of each species.
        {
            energy[s] = energy[s] * forwardReactionTemperature + thermodynamicProperties[s].hs0;
        }
    }
    forwardEnergy = density * mixGasCv * forwardReactionTemperature;
    backwardEnergy = density * mixGasCv * backwardReactionTemperature;

    for (int r = 0; r < numberOfReaction; ++ r)
    {
        forwardRates[r] = 1.0;
        backwardRates[r] = 1.0;
        collisionRates[r] = 0.0;
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            densityPerMass = densityMassFraction[s] * speciesMassReciprocal[s];
            for (int i = 0; i < forwardStoichNumber[r][s]; ++ i)
            {
                forwardRates[r] *= densityPerMass;
            }
            for (int i = 0; i < backwardStoichNumber[r][s]; ++ i)
            {
                backwardRates[r] *= densityPerMass;
            }
            collisionRates[r] += densityPerMass * thirdBodyStoichNumber[r][s];
        }
        if (fabs(collisionRates[r]) < SMALL) // RCr
        {
            collisionRates[r] = 1.0;
        }

        //! Compute equilibrium constant of the front reaction.
        kf = forwardRateCoefA[r] * pow(forwardReactionTemperature, forwardRateCoefB[r]) * exp(-forwardRateCoefC[r] / forwardReactionTemperature);
        //! Compute equilibrium constant of the reverse reaction.
        kb = backwardRateCoefA[r] * pow(backwardReactionTemperature, backwardRateCoefB[r]) * exp(-backwardRateCoefC[r] / backwardReactionTemperature);
        forwardRates[r] = kf * forwardRates[r];     //! RAr.
        backwardRates[r] = kb * backwardRates[r];   //! RBr.

        Drf = (forwardRateCoefB[r] + forwardRateCoefC[r] / forwardReactionTemperature) / forwardEnergy;
        Drb = (backwardRateCoefB[r] + backwardRateCoefC[r] / backwardReactionTemperature) / backwardEnergy;

        //! Compute the partial derivatives to each species, (dRr / dCs).
        for (int s = 0; s < this->nLastSpeciesIndex; ++ s)
        {
            //! One-Temperature Model.
            coefA = ionType[s] * speciesMassReciprocal[s] * this->electronMassDimension;
            gamaj = energy[nLastSpeciesIndex] - energy[s] + coefA * (energy[nLastSpeciesIndex] - energy[electronIndex]);

            theta[0]  = forwardStoichNumber[r][s] / densityMassFraction[s] - forwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex];
            theta[0] += coefA * (forwardStoichNumber[r][electronIndex] / densityMassFraction[electronIndex] - forwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex]);
            theta[0] += Drf * gamaj;
            theta[0]  = forwardRates[r] * theta[0];

            theta[1]  = backwardStoichNumber[r][s] / densityMassFraction[s] - backwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex];
            theta[1] += coefA * (backwardStoichNumber[r][electronIndex] / densityMassFraction[electronIndex] - backwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex]);
            theta[1] += Drb * gamaj;
            theta[1]  = backwardRates[r] * theta[1];

            theta[2]  = thirdBodyStoichNumber[r][s] * speciesMassReciprocal[s] - thirdBodyStoichNumber[r][nLastSpeciesIndex] * speciesMassReciprocal[nLastSpeciesIndex];
            theta[2] += coefA * (thirdBodyStoichNumber[r][electronIndex] * speciesMassReciprocal[electronIndex] - thirdBodyStoichNumber[r][nLastSpeciesIndex] * speciesMassReciprocal[nLastSpeciesIndex]);

            //! Compute the partial derivatives to the species j. (dRr / dCs).
            partialValues[r][nm + s] = (theta[0] - theta[1]) * collisionRates[r] + (forwardRates[r] - backwardRates[r]) * theta[2];
        }

        //! Compute the partial derivative to vibration and electron energy.
        for (int it = 1; it < ntmodel; ++ it)
        {
            phir = (Drf * forwardRates[r] - Drb * backwardRates[r]) * collisionRates[r]; // (dRr / dEv).
            partialValues[r][nl + it] = -phir;
        }

    }

    //! Compute the chemical source term.
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        chemicalSourceTerms[s] = 0.0;
        for (int r = 0; r < numberOfReaction; ++ r)
        {
            chemicalSourceTerms[s] += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * (forwardRates[r] - backwardRates[r]) * collisionRates[r];
        }
        chemicalSourceTerms[s] *= referenceSpeciesEquation * speciesMass[s] * cellVolume;
    }

    //! Construct elements of the Jacobian matrix of the chemical source term.
    for (int s = 0; s < this->nLastSpeciesIndex; ++ s)
    {
        sourceDerivatives[s][s] = 0.0;
        for (int r = 0; r < numberOfReaction; ++ r)
        {
            sourceDerivatives[s][s] += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * partialValues[r][nm + s]; // (dRr / dcs).
        }
        sourceDerivatives[s][s] *= reference1 * speciesMass[s] * cellVolume;// Translate to the nondimensional values.
    }

    for (int it = 1; it < ntmodel; ++ it)      //! Compute the partial derivatives of species source term to the vibration energy.
    {
        for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
        {
            sourceDerivatives[s][numberOfSpecies + it - 1] = 0.0;
            for (int r = 0; r < numberOfReaction; ++ r)
            {
                sourceDerivatives[s][numberOfSpecies + it - 1] += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * partialValues[r][nl + it]; // (dRr / dEv).
            }
            sourceDerivatives[s][numberOfSpecies + it - 1] *= reference2 * speciesMass[s] * cellVolume;
        }
    }

    delete []energy;
    delete []forwardRates;
    delete []backwardRates;
    delete []collisionRates;
    delete []densityMassFraction;
    DelPointer2(partialValues);
}

void MixedGas::ComputeChemicalSourceAndJacobian(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{
    if (nIsChemicalFreeze > 0)
    {
        for (int iCell = 0; iCell < nCellNumber; ++ iCell)
        {
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                chemicalSourceTerms[s][iCell] = 0.0;
                sourceDerivatives[s][s][iCell] = 0.0;
            }
        }
        return;
    }
    using namespace IDX;     //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();

    //! For source term in species transport equation.
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble reference1 = referenceDensity * referenceSpeciesEquation;
    RDouble reference2 = referenceVelocity * reynoldsReferenceLengthDimensional;

    RDouble *forwardRateCoefA = reactionRate->GetForwardReactionRateCoefficientA();
    RDouble *forwardRateCoefB = reactionRate->GetForwardReactionRateCoefficientB();
    RDouble *forwardRateCoefC = reactionRate->GetForwardReactionRateCoefficientC();
    RDouble *backwardRateCoefA = reactionRate->GetBackwardReactionRateCoefficientA();
    RDouble *backwardRateCoefB = reactionRate->GetBackwardReactionRateCoefficientB();
    RDouble *backwardRateCoefC = reactionRate->GetBackwardReactionRateCoefficientC();

    int *reactionType = stoichiometricEquation->GetChemicalReactionType();
    bool *isCollisionFlag = stoichiometricEquation->GetCollisionReactionFlag();

    int **forwardStoichNumber = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
    int **backwardStoichNumber = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
    RDouble **thirdBodyStoichNumber = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();

    //! Record the coefficients in the chemical reactions.
    RDouble *forwardRates = new RDouble [numberOfReaction];
    RDouble *backwardRates = new RDouble [numberOfReaction];
    RDouble *collisionRates = new RDouble [numberOfReaction];
    RDouble *specificHeatCv = new RDouble [numberOfSpecies];
    RDouble *energy = new RDouble [numberOfSpecies];
    RDouble *densityMassFraction = new RDouble [numberOfSpecies];
    int nEquation = this->nm + this->numberOfSpecies + this->ntmodel - 1;
    RDouble **partialValues = NewPointer2 <RDouble> (numberOfReaction, nEquation);
    //RDouble primVars[MAX_SPECIES_NUM], speciesCvtr[MAX_SPECIES_NUM], speciesCvv[MAX_SPECIES_NUM], totalCvtr, totalCvv;
    //RDouble speciesEtr[MAX_SPECIES_NUM], speciesEv[MAX_SPECIES_NUM], speciesH[MAX_SPECIES_NUM], totalEtr, totalEv, totalH0;
    //The power of the translational temperature in Park model.
    RDouble Tc1[6] = {parkVDPower, 0.0, 1.0, 1.0, 1.0, 0.0}; //for the forward reaction.
    RDouble Tc2[6] = {1.0, 0.5, 1.0, 1.0, 0.5, 0.0}; //for the backward reaction.

    //! The dimensional value of control temperature of chemical reaction;
    RDouble forwardReactionTemperature, backwardReactionTemperature, Tf[6], Tb[6], trsTemperature, vibTemperature;
    RDouble densityPerMass, density, mixGasCv, kf, kb, forwardEnergy, backwardEnergy;
    RDouble gamaj = 0.0, Drf = 0.0, Drb = 0.0, phir = 0.0, theta[3] = {0.0}, coefA = 0.0;
    int electronIndex = numberOfSpecies - 1; // The index of electron in the array of species.
    if (nElectronIndex >= 0)
    {
        electronIndex = nElectronIndex;
    }

    //! Compute the elements of Jacobian matrix.
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();

    for (int iCell = 0; iCell < nCellNumber; ++ iCell)
    {
        //Obtain the dimensional temperature.
        trsTemperature = MAX(temperatures[ITT][iCell] * referenceTemperature, SMALL);
        vibTemperature = trsTemperature;
        if (ntmodel > 1)
        {
            vibTemperature = MAX(temperatures[ITV][iCell] * referenceTemperature, SMALL);
        }
        //the control temperature of dissociation reaction, e.g. AB = A + B.
        Tf[0] = pow(trsTemperature, parkVDPower) * pow(vibTemperature, 1.0 - parkVDPower);
        Tb[0] = trsTemperature;
        //the control temperature of electron collision dissociation reaction, e.g. AB + e- = A + B + e-.
        Tf[1] = vibTemperature;
        Tb[1] = pow(trsTemperature, 0.5) * pow(vibTemperature, 0.5);
        //the control temperature of neutral exchange reaction, e.g. AB + C = A + BC.
        Tf[2] = trsTemperature;
        Tb[2] = trsTemperature;
        //the control temperature of charge exchange reaction, e.g. AB+ + C = A+ + BC.
        Tf[3] = trsTemperature;
        Tb[3] = trsTemperature;
        //the control temperature of recombination ionization reaction, e.g. A + B = AB+ + e-.
        Tf[4] = trsTemperature;
        Tb[4] = Tb[1];
        //the control temperature of electron collision ionization reaction, e.g. A + e- = A+ + 2e-.
        Tf[5] = vibTemperature;
        Tb[5] = vibTemperature;

        //! The control temperature of front reaction.
        forwardReactionTemperature = Tf[0];
        //! The control temperature of back reaction. it is temperature of translation and rotation.
        backwardReactionTemperature = Tb[0];
        density = primitiveVariables[IR][iCell] * referenceDensity;
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            massFractions[s] = primitiveVariables[nm + s][iCell];
            densityMassFraction[s] = MAX(SMALL, density * massFractions[s]);
        }

        if (ntmodel == 1)      //! One-Temperature Model.
        {
            //! Compute the specific heat at constant volume of each species.
            this->ComputeSpeciesConstantVolumeSpecificHeatDimensional(temperatures[ITT][iCell], specificHeatCv);
            //this->ComputeSpeciesEnergy(temperatures[ITT][iCell], energy);
            this->ComputeMixtureByMassFraction(massFractions, specificHeatCv, mixGasCv);
            for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute internal energy of each species.
            {
                energy[s] = specificHeatCv[s] * forwardReactionTemperature;
            }
        }
        else         //! Multi-Temperature Model.
        {
            //if (ntmodel == 2 && nTEnergyModel == 1)
            //{
            //    for (int s = 0; s < numberOfSpecies; ++ s)
            //    {
            //        primVars[nm + s] = massFractions[s];
            //    }
            //    RDouble tArray[3] = {temperatures[ITT][iCell], temperatures[ITV][iCell], temperatures[ITV][iCell]};
            //    this->ComputeSpecificHeat(primVars, tArray, speciesCvtr, speciesCvv, totalCvtr, totalCvv);
            //    totalCvtr = totalCvtr / referenceSpecificHeat;
            //    totalCvv  = totalCvv / referenceSpecificHeat;

            //    this->ComputeEnergy(primVars, tArray, speciesEtr, speciesEv, speciesH, totalEtr, totalEv, totalH0);
            //    for (int s = 0; s < numberOfSpecies; ++ s)
            //    {
            //        energy[s] = speciesEtr[s] * referenceEnergy;
            //    }
            //}
            //else
            //{
                // Compute the specific heat at constant volume of each species.
                this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(specificHeatCv);
                this->ComputeMixtureByMassFraction(massFractions, specificHeatCv, mixGasCv);
                for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute internal energy of each species.
                {
                    energy[s] = specificHeatCv[s] * forwardReactionTemperature + thermodynamicProperties[s].hs0;
                    //energy[s] = specificHeatCv[s] * trsTemperature + thermodynamicProperties[s].hs0; //this is correct.
                }
            //}
        }

        forwardEnergy = density * mixGasCv * forwardReactionTemperature;
        backwardEnergy = density * mixGasCv * backwardReactionTemperature;

        for (int r = 0; r < numberOfReaction; ++ r)
        {
            forwardRates[r] = 1.0;
            backwardRates[r] = 1.0;
            collisionRates[r] = 0.0;
            int nType = reactionType[r]; //The type of chemical reaction.
            forwardReactionTemperature = Tf[nType];
            backwardReactionTemperature = Tb[nType];
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                densityPerMass = density * massFractions[s] * speciesMassReciprocal[s];
                for (int i = 0; i < forwardStoichNumber[r][s]; ++ i)
                {
                    forwardRates[r] *= densityPerMass;
                }
                for (int i = 0; i < backwardStoichNumber[r][s]; ++ i)
                {
                    backwardRates[r] *= densityPerMass;
                }
                collisionRates[r] += densityPerMass * thirdBodyStoichNumber[r][s];
            }
            if (!isCollisionFlag[r]) // RCr
            {
                collisionRates[r] = 1.0;
            }

            //! Compute equilibrium constant of the front reaction.
            kf =forwardRateCoefA[r] * pow(forwardReactionTemperature, forwardRateCoefB[r]) * exp(-forwardRateCoefC[r] / forwardReactionTemperature);
            //! Compute equilibrium constant of the reverse reaction.
            kb =backwardRateCoefA[r] * pow(backwardReactionTemperature, backwardRateCoefB[r]) * exp(-backwardRateCoefC[r] / backwardReactionTemperature);
            forwardRates[r] = kf * forwardRates[r];     //! RAr.
            backwardRates[r] = kb * backwardRates[r];   //! RBr.

            //Compute the partial derivative coefficient of the reaction temperature to the vibrational energy, namely dT / d(rho * Ev);
            //if (ntmodel == 2 && nTEnergyModel == 1)
            //{
            //    RDouble TvDivT = vibTemperature / trsTemperature;
            //    RDouble phi_1 = Tc1[nType] * pow(TvDivT, 1.0 - Tc1[nType]) / (density * totalCvtr);
            //    phi_1 -= (1.0 - Tc1[nType]) * pow (TvDivT, -Tc1[nType]) / (density * totalCvv);
            //    Drf = phi_1 * (forwardRateCoefB[r] + forwardRateCoefC[r] / forwardReactionTemperature) / forwardReactionTemperature;

            //    RDouble phi_2 = Tc2[nType] * pow(TvDivT, 1.0 - Tc2[nType]) / (density * totalCvtr);
            //    phi_2 -= (1.0 - Tc2[nType]) * pow (TvDivT, -Tc2[nType]) / (density * totalCvv);
            //    Drb = phi_2 * (backwardRateCoefB[r] + backwardRateCoefC[r] / backwardReactionTemperature) / backwardReactionTemperature;
            //}
            //else
            //{
                Drf = (forwardRateCoefB[r] + forwardRateCoefC[r] / forwardReactionTemperature) / forwardEnergy;
                Drb = (backwardRateCoefB[r] + backwardRateCoefC[r] / backwardReactionTemperature) / backwardEnergy;
            //}

            //! Compute the partial derivatives to each species, (dRr / dCs).
            for (int s = 0; s < this->nLastSpeciesIndex; ++ s)
            {
                //! One-Temperature Model.
                coefA = ionType[s] * speciesMassReciprocal[s] * this->electronMassDimension;
                gamaj = energy[nLastSpeciesIndex] - energy[s] + coefA * (energy[nLastSpeciesIndex] - energy[electronIndex]);

                theta[0]  = forwardStoichNumber[r][s] / densityMassFraction[s] - forwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex];
                theta[0] += coefA * (forwardStoichNumber[r][electronIndex] / densityMassFraction[electronIndex] - forwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex]);
                theta[0] += Drf * gamaj;
                theta[0]  = forwardRates[r] * theta[0];

                theta[1]  = backwardStoichNumber[r][s] / densityMassFraction[s] - backwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex];
                theta[1] += coefA * (backwardStoichNumber[r][electronIndex] / densityMassFraction[electronIndex] - backwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex]);
                theta[1] += Drb * gamaj;
                theta[1]  = backwardRates[r] * theta[1];

                theta[2]  = thirdBodyStoichNumber[r][s] * speciesMassReciprocal[s] - thirdBodyStoichNumber[r][nLastSpeciesIndex] * speciesMassReciprocal[nLastSpeciesIndex];
                theta[2] += coefA * (thirdBodyStoichNumber[r][electronIndex] * speciesMassReciprocal[electronIndex] - thirdBodyStoichNumber[r][nLastSpeciesIndex] * speciesMassReciprocal[nLastSpeciesIndex]);

                //! Compute the partial derivatives to the species j. (dRr / dCs).
                partialValues[r][nm + s] = (theta[0] - theta[1]) * collisionRates[r] + (forwardRates[r] - backwardRates[r]) * theta[2];
            }

            //! Compute the partial derivative to vibration and electron energy.
            for (int it = 1; it < ntmodel; ++ it)
            {
                phir = (Drf * forwardRates[r] - Drb * backwardRates[r]) * collisionRates[r]; // (dRr / dEv).
                partialValues[r][nl + it] = -phir;
            }

        }

        RDouble fSum = 0.0;
        //! Compute the chemical source term.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            fSum = 0.0;
            for (int r = 0; r < numberOfReaction; ++ r)
            {
                fSum += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * (forwardRates[r] - backwardRates[r]) * collisionRates[r];
            }
            chemicalSourceTerms[s][iCell] = referenceSpeciesEquation * speciesMass[s] * fSum * cellVolume[iCell];
        }

        //! Construct elements of the Jacobian matrix of the chemical source term.
        for (int s = 0; s < this->nLastSpeciesIndex; ++ s)
        {
            fSum = 0.0;
            for (int r = 0; r < numberOfReaction; ++ r)
            {
                fSum += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * partialValues[r][nm + s]; // (dRr / dcs).
            }
            sourceDerivatives[s][s][iCell] = reference1 * speciesMass[s] * fSum * cellVolume[iCell];// Translate to the nondimensional values.  //?
        }

        for (int it = 1; it < ntmodel; ++ it)      //! Compute the partial derivatives of species source term to the vibration energy.
        {
            for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
            {
                fSum = 0.0;
                for (int r = 0; r < numberOfReaction; ++ r)
                {
                    fSum += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * partialValues[r][nl + it]; // (dRr / dEv).
                }
                sourceDerivatives[s][numberOfSpecies + it - 1][iCell] = reference2 * speciesMass[s] * fSum * cellVolume[iCell];               //?
            }
        }

        //Inside of the loop.
    }

    delete []energy;
    delete []specificHeatCv;
    delete []forwardRates;
    delete []backwardRates;
    delete []collisionRates;
    delete []densityMassFraction;
    DelPointer2(partialValues);
}


void MixedGas::ComputeChemicalSourceAndJacobianD3(int *Inzone, RDouble4D &primitiveVariables, RDouble4D &temperatures,RDouble4D &speciesEs, RDouble4D &speciesCvs, RDouble3D &totalCv, RDouble3D &cellVolume, RDouble4D &chemicalSourceTerms, RDouble5D &sourceDerivatives)
{
    using namespace IDX;     //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();

    //! For source term in species transport equation.
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble reference1 = referenceDensity * referenceSpeciesEquation;
    RDouble reference2 = referenceVelocity * reynoldsReferenceLengthDimensional;

    RDouble *forwardRateCoefA = reactionRate->GetForwardReactionRateCoefficientA();
    RDouble *forwardRateCoefB = reactionRate->GetForwardReactionRateCoefficientB();
    RDouble *forwardRateCoefC = reactionRate->GetForwardReactionRateCoefficientC();
    RDouble *backwardRateCoefA = reactionRate->GetBackwardReactionRateCoefficientA();
    RDouble *backwardRateCoefB = reactionRate->GetBackwardReactionRateCoefficientB();
    RDouble *backwardRateCoefC = reactionRate->GetBackwardReactionRateCoefficientC();

    int *reactionType = stoichiometricEquation->GetChemicalReactionType();
    bool *isCollisionFlag = stoichiometricEquation->GetCollisionReactionFlag();

    int **forwardStoichNumber = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
    int **backwardStoichNumber = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
    RDouble **thirdBodyStoichNumber = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();

    int nEquation = this->nm + this->numberOfSpecies + this->ntmodel - 1;

    //! Record the coefficients in the chemical reactions.
    RDouble *forwardRates = Thermo_Energy_temparay.forwardRates;
    RDouble *backwardRates = Thermo_Energy_temparay.backwardRates;
    RDouble *collisionRates = Thermo_Energy_temparay.collisionRates;
    RDouble *specificHeatCv = Thermo_Energy_temparay.Cvs;
    RDouble *energy = Thermo_Energy_temparay.Es;
    RDouble *densityMassFraction = Thermo_Energy_temparay.Cvvs;
    RDouble **partialValues =  Thermo_Energy_temparay.partialValues;
    RDouble **partialValues2 =  Thermo_Energy_temparay.partialValues2;
    //RDouble primVars[MAX_SPECIES_NUM], speciesCvtr[MAX_SPECIES_NUM], speciesCvv[MAX_SPECIES_NUM], totalCvtr, totalCvv;
    //RDouble speciesEtr[MAX_SPECIES_NUM], speciesEv[MAX_SPECIES_NUM], speciesH[MAX_SPECIES_NUM], totalEtr, totalEv, totalH0;
    //The power of the translational temperature in Park model.
    RDouble Tc1[6] = {parkVDPower, 0.0, 1.0, 1.0, 1.0, 0.0}; //for the forward reaction.
    RDouble Tc2[6] = {1.0, 0.5, 1.0, 1.0, 0.5, 0.0}; //for the backward reaction.

    //! The dimensional value of control temperature of chemical reaction;
    RDouble forwardReactionTemperature, backwardReactionTemperature, Tf[6], Tb[6], trsTemperature, vibTemperature;
    RDouble densityPerMass, density, mixGasCv, kf, kb, forwardEnergy, backwardEnergy;
    RDouble gamaj = 0.0, Drf = 0.0, Drb = 0.0, phir = 0.0, theta[3] = {0.0}, coefA = 0.0;
    int electronIndex = numberOfSpecies - 1; // The index of electron in the array of species.
    if (nElectronIndex >= 0)
    {
        electronIndex = nElectronIndex;
    }

    //! Compute the elements of Jacobian matrix.
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();

    int ist, ied, jst, jed, kst, ked;

    ist = Inzone[0];
    ied = Inzone[1];
    jst = Inzone[2];
    jed = Inzone[3];
    kst = Inzone[4];
    ked = Inzone[5];

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                trsTemperature = MAX(temperatures(i, j, k, ITT) * referenceTemperature, SMALL);
                vibTemperature = trsTemperature;
                if (ntmodel > 1)
                {
                    vibTemperature = MAX(temperatures(i, j, k, ITV)  * referenceTemperature, SMALL);
                }
                //the control temperature of dissociation reaction, e.g. AB = A + B.
                Tf[0] = pow(static_cast<RealData>(trsTemperature), static_cast<RealData>(parkVDPower)) * pow(static_cast<RealData>(vibTemperature), static_cast<RealData>(1.0 - parkVDPower));
                Tb[0] = trsTemperature;
                //the control temperature of electron collision dissociation reaction, e.g. AB + e- = A + B + e-.
                Tf[1] = vibTemperature;
                Tb[1] = sqrt(static_cast<RealData>(trsTemperature * vibTemperature));
                //the control temperature of neutral exchange reaction, e.g. AB + C = A + BC.
                Tf[2] = trsTemperature;
                Tb[2] = trsTemperature;
                //the control temperature of charge exchange reaction, e.g. AB+ + C = A+ + BC.
                Tf[3] = trsTemperature;
                Tb[3] = trsTemperature;
                //the control temperature of recombination ionization reaction, e.g. A + B = AB+ + e-.
                Tf[4] = trsTemperature;
                Tb[4] = Tb[1];
                //the control temperature of electron collision ionization reaction, e.g. A + e- = A+ + 2e-.
                Tf[5] = vibTemperature;
                Tb[5] = vibTemperature;

                //! The control temperature of front reaction.
                forwardReactionTemperature = Tf[0];
                //! The control temperature of back reaction. it is temperature of translation and rotation.
                backwardReactionTemperature = Tb[0];
                density = primitiveVariables(i, j, k, IR) * referenceDensity;
                for (int s = 0; s < numberOfSpecies; ++ s)
                {
                    massFractions[s] = primitiveVariables(i, j, k, nm + s);
                    densityMassFraction[s] = MAX(SMALL, density * massFractions[s]);

                    specificHeatCv[s] = speciesCvs(i, j, k, s) / referenceSpecificHeat;
                }
                mixGasCv = totalCv(i, j, k) / referenceSpecificHeat;


                if(nChemcalSourceEsMethod==0)
                {
                    if (ntmodel == 1)      //! One-Temperature Model.
                    {
                        for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute internal energy of each species.
                        {
                            energy[s] = specificHeatCv[s] * forwardReactionTemperature;
                        }
                    }
                    else         //! Multi-Temperature Model.
                    {
                        for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute internal energy of each species.
                        {
                            energy[s] = specificHeatCv[s] * forwardReactionTemperature + thermodynamicProperties[s].hs0;
                        }
                    }
                }
                else
                {
                    for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute internal energy of each species.
                    {
                            energy[s] = speciesEs(i, j, k, s) * referenceEnergy;
                    }
                
                }
                

                forwardEnergy = density * mixGasCv * forwardReactionTemperature;
                backwardEnergy = density * mixGasCv * backwardReactionTemperature;

                for (int r = 0; r < numberOfReaction; ++ r)
                {
                    forwardRates[r] = 1.0;
                    backwardRates[r] = 1.0;
                    collisionRates[r] = 0.0;
                    int nType = reactionType[r]; //The type of chemical reaction.
                    forwardReactionTemperature = Tf[nType];
                    backwardReactionTemperature = Tb[nType];
                    for (int s = 0; s < numberOfSpecies; ++ s)
                    {
                        densityPerMass = density * massFractions[s] * speciesMassReciprocal[s];
                        for (int i = 0; i < forwardStoichNumber[r][s]; ++ i)
                        {
                            forwardRates[r] *= densityPerMass;
                        }
                        for (int i = 0; i < backwardStoichNumber[r][s]; ++ i)
                        {
                            backwardRates[r] *= densityPerMass;
                        }
                        collisionRates[r] += densityPerMass * thirdBodyStoichNumber[r][s];
                    }
                    if (!isCollisionFlag[r]) // RCr
                    {
                        collisionRates[r] = 1.0;
                    }

                    //! Compute equilibrium constant of the front reaction.
                    //kf =forwardRateCoefA[r]  * pow(static_cast<RealData>(forwardReactionTemperature),  static_cast<RealData>(forwardRateCoefB[r]))  * exp(static_cast<RealData>(-forwardRateCoefC[r] / forwardReactionTemperature));
                    //! Compute equilibrium constant of the reverse reaction.
                    //kb =backwardRateCoefA[r] * pow(static_cast<RealData>(backwardReactionTemperature), static_cast<RealData>(backwardRateCoefB[r])) * exp(static_cast<RealData>(-backwardRateCoefC[r] / backwardReactionTemperature));
                    kf = min(lnMaximum,log(static_cast<RealData>(forwardRateCoefA[r]))  +   forwardRateCoefB[r]  * log(static_cast<RealData>(forwardReactionTemperature))  - forwardRateCoefC[r]  / forwardReactionTemperature);
                    kb = min(lnMaximum,log(static_cast<RealData>(backwardRateCoefA[r])) +   backwardRateCoefB[r] * log(static_cast<RealData>(backwardReactionTemperature)) - backwardRateCoefC[r] / backwardReactionTemperature);

                    kf = exp(static_cast<RealData>(kf));
                    kb = exp(static_cast<RealData>(kb));
                    
                   /* if(r>10)
                    {
                        kf=0.1*kf;
                        kb=0.1*kb;                    
                    }*/

                    forwardRates[r] = kf * forwardRates[r];     //! RAr.
                    backwardRates[r] = kb * backwardRates[r];   //! RBr.

                    //Compute the partial derivative coefficient of the reaction temperature to the vibrational energy, namely dT / d(rho * Ev);
                    //if (ntmodel == 2 && nTEnergyModel == 1)
                    //{
                    //    RDouble TvDivT = vibTemperature / trsTemperature;
                    //    RDouble phi_1 = Tc1[nType] * pow(TvDivT, 1.0 - Tc1[nType]) / (density * totalCvtr);
                    //    phi_1 -= (1.0 - Tc1[nType]) * pow (TvDivT, -Tc1[nType]) / (density * totalCvv);
                    //    Drf = phi_1 * (forwardRateCoefB[r] + forwardRateCoefC[r] / forwardReactionTemperature) / forwardReactionTemperature;

                    //    RDouble phi_2 = Tc2[nType] * pow(TvDivT, 1.0 - Tc2[nType]) / (density * totalCvtr);
                    //    phi_2 -= (1.0 - Tc2[nType]) * pow (TvDivT, -Tc2[nType]) / (density * totalCvv);
                    //    Drb = phi_2 * (backwardRateCoefB[r] + backwardRateCoefC[r] / backwardReactionTemperature) / backwardReactionTemperature;
                    //}
                    //else
                    //{
                        Drf = (forwardRateCoefB[r] + forwardRateCoefC[r] / forwardReactionTemperature) / forwardEnergy;
                        Drb = (backwardRateCoefB[r] + backwardRateCoefC[r] / backwardReactionTemperature) / backwardEnergy;
                    //}

                    //! Compute the partial derivatives to each species, (dRr / dCs).
                    for (int s = 0; s < this->nLastSpeciesIndex; ++ s)
                    {
                        //! One-Temperature Model.
                        coefA = ionType[s] * speciesMassReciprocal[s] * this->electronMassDimension;
                        gamaj = energy[nLastSpeciesIndex] - energy[s] + coefA * (energy[nLastSpeciesIndex] - energy[electronIndex]);

                        theta[0]  = forwardStoichNumber[r][s] / densityMassFraction[s] - forwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex];
                        theta[0] += coefA * (forwardStoichNumber[r][electronIndex] / densityMassFraction[electronIndex] - forwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex]);
                        //theta[0] += Drf * gamaj;
                        theta[0]  = forwardRates[r] * theta[0];

                        theta[1]  = backwardStoichNumber[r][s] / densityMassFraction[s] - backwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex];
                        theta[1] += coefA * (backwardStoichNumber[r][electronIndex] / densityMassFraction[electronIndex] - backwardStoichNumber[r][nLastSpeciesIndex] / densityMassFraction[nLastSpeciesIndex]);
                        //theta[1] += Drb * gamaj;
                        theta[1]  = backwardRates[r] * theta[1];

                        theta[2]  = thirdBodyStoichNumber[r][s] * speciesMassReciprocal[s] - thirdBodyStoichNumber[r][nLastSpeciesIndex] * speciesMassReciprocal[nLastSpeciesIndex];
                        theta[2] += coefA * (thirdBodyStoichNumber[r][electronIndex] * speciesMassReciprocal[electronIndex] - thirdBodyStoichNumber[r][nLastSpeciesIndex] * speciesMassReciprocal[nLastSpeciesIndex]);

                        //! Compute the partial derivatives to the species j. (dRr / dCs).
                        partialValues[r][nm + s] = (theta[0] - theta[1]) * collisionRates[r] + (forwardRates[r] - backwardRates[r]) * theta[2];
                        partialValues2[r][nm + s] = (forwardRates[r] * Drf - backwardRates[r] * Drb)* gamaj * collisionRates[r];
                    }

                    //! Compute the partial derivative to vibration and electron energy.
                    for (int it = 1; it < ntmodel; ++ it)
                    {
                        phir = (Drf * forwardRates[r] - Drb * backwardRates[r]) * collisionRates[r]; // (dRr / dEv).
                        partialValues[r][nl + it] = -phir;
                    }

                }

                RDouble fSum = 0.0, fSum1 = 0.0;
                //! Compute the chemical source term.
                for (int s = 0; s < numberOfSpecies; ++ s)
                {
                    fSum = 0.0;
                    for (int r = 0; r < numberOfReaction; ++ r)
                    {
                        fSum += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * (forwardRates[r] - backwardRates[r]) * collisionRates[r];
                    }
                    chemicalSourceTerms(i, j, k, s) = referenceSpeciesEquation * speciesMass[s] * fSum * cellVolume(i, j, k);
                }

                //! Construct elements of the Jacobian matrix of the chemical source term.
                for (int s = 0; s < this->nLastSpeciesIndex; ++ s)
                {
                    fSum = 0.0;
                    for (int r = 0; r < numberOfReaction; ++ r)
                    {
                        fSum += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * partialValues[r][nm + s]; // (dRr / dcs).
                    }

                    fSum1 = 0.0;
                    for (int r = 0; r < numberOfReaction; ++ r)
                    {
                        fSum1 += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * partialValues2[r][nm + s]; // (dRr / dcs).
                    }
                   //
                    if(nChemcalSourceModified==0)
                    {
                        fSum = MIN(0.0, fSum + fSum1);
                    }
                    else if(nChemcalSourceModified==1)
                    {
                        fSum = MIN(0.0, fSum) + MIN(0.0, fSum1);
                    }
                    else
                    {
                        fSum = -abs(fSum) - abs(fSum1);
                    }
                    sourceDerivatives(i, j, k, s, s) = reference1 * speciesMass[s] * fSum * cellVolume(i, j, k);// Translate to the nondimensional values.  //?
                }

                for (int it = 1; it < ntmodel; ++ it)      //! Compute the partial derivatives of species source term to the vibration energy.
                {
                    for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
                    {
                        fSum = 0.0;
                        for (int r = 0; r < numberOfReaction; ++ r)
                        {
                            fSum += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * partialValues[r][nl + it]; // (dRr / dEv).
                        }
                        sourceDerivatives(i, j, k, s, numberOfSpecies + it - 1) = reference2 * speciesMass[s] * fSum * cellVolume(i, j, k);               //?
                    }
                }

                /*double r0[10],r1[10];

                if(i==114&&j==29&&k==1)
                {
                    for (int s = 0; s < numberOfSpecies; ++ s)
                    {
                        r0[s] = chemicalSourceTerms(i, j, k, s);
                    }


                    for (int s = 0; s < numberOfSpecies; ++ s)
                    {
                        r1[s] = sourceDerivatives(i, j, k, s, s);
                    }

                    r1[0] = 0;
                }*/


        //Inside of the loop.
            }
        }
    }

    //delete []energy;
    //delete []specificHeatCv;
    //delete []forwardRates;
    //delete []backwardRates;
    //delete []collisionRates;
    //delete []densityMassFraction;
    //DelPointer2(partialValues);
}



void MixedGas::ComputeChemicalSourceDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{
    //! The function is used for chemical reaction without electron.

    using namespace IDX;     //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();

    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();

    //! The dimensional value of control temperature of chemical reaction.
    RDouble forwardReactionTemperature, backwardReactionTemperature;
    RDouble densityPerMass, density, mixGasCv, squareVelocity;

    //! For source term in species transport equation.
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble reference1 = referenceDensity * referenceSpeciesEquation;

    RDouble *forwardRateCoefA = reactionRate->GetForwardReactionRateCoefficientA();
    RDouble *forwardRateCoefB = reactionRate->GetForwardReactionRateCoefficientB();
    RDouble *forwardRateCoefC = reactionRate->GetForwardReactionRateCoefficientC();

    RDouble *backwardRateCoefA = reactionRate->GetBackwardReactionRateCoefficientA();
    RDouble *backwardRateCoefB = reactionRate->GetBackwardReactionRateCoefficientB();
    RDouble *backwardRateCoefC = reactionRate->GetBackwardReactionRateCoefficientC();

    int **forwardStoichNumber = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
    int **backwardStoichNumber = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
    RDouble **thirdBodyStoichNumber = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();

    //! Record the coefficients in the chemical reactions.
    RDouble *forwardRates = new RDouble [numberOfReaction];
    RDouble *backwardRates = new RDouble [numberOfReaction];
    RDouble *collisionRates = new RDouble [numberOfReaction];
    RDouble *energy = new RDouble [numberOfSpecies];
    RDouble *massFraction = new RDouble [numberOfSpecies];
    RDouble *densityMassFraction = new RDouble [numberOfSpecies];
    int nEquation = this->nm + this->numberOfSpecies + this->ntmodel - 1;
    RDouble **partialValues = NewPointer2 <RDouble> (numberOfReaction, nEquation);

    //! Compute the elements of Jacobian matrix.
    RDouble alpha = 0.0, gamaj = 0.0, Drf = 0.0, Drb = 0.0, phir = 0.0, theta[3] = {0.0};
    RDouble uVelocity, vVelocity, wVelocity;
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();
    for (int iCell = 0; iCell < nCellNumber; ++ iCell)
    {
        //! The control temperature of front reaction.
        forwardReactionTemperature = MAX(temperatures[ITT][iCell] * referenceTemperature, SMALL);
        //! The control temperature of back reaction. it is temperature of translation and rotation.
        backwardReactionTemperature = forwardReactionTemperature;

        density = primitiveVariables[IR][iCell] * referenceDensity;
        uVelocity = primitiveVariables[IU][iCell] * referenceVelocity;
        vVelocity = primitiveVariables[IV][iCell] * referenceVelocity;
        wVelocity = primitiveVariables[IW][iCell] * referenceVelocity;
        squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            massFraction[s] = primitiveVariables[nm + s][iCell];
            densityMassFraction[s] = MAX(SMALL, density * massFraction[s]);
        }

        if (ntmodel == 1)     //! One-Temperature Model.
        {
            //! Compute the specific heat at constant volume of each species.
            this->ComputeSpeciesConstantVolumeSpecificHeatDimensional(temperatures[ITT][iCell], energy);
            this->ComputeMixtureByMassFraction(massFraction, energy, mixGasCv);
            for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute internal energy of each species.
            {
                energy[s] *= forwardReactionTemperature;
            }
        }
        else     //! Multi-Temperature Model.
        {
            //! Compute the specific heat at constant volume of each species.
            this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(energy);
            this->ComputeMixtureByMassFraction(massFraction, energy, mixGasCv);
            for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute internal energy of each species.
            {
                energy[s] = energy[s] * forwardReactionTemperature + thermodynamicProperties[s].hs0;
            }
        }

        for (int r = 0; r < numberOfReaction; ++ r)
        {
            //! Compute equilibrium constant of the front reaction.
            RDouble kf = forwardRateCoefA[r] * pow(forwardReactionTemperature, forwardRateCoefB[r]) * exp(-forwardRateCoefC[r] / forwardReactionTemperature);
            //! Compute equilibrium constant of the reverse reaction.
            RDouble kb = backwardRateCoefA[r] * pow(backwardReactionTemperature, backwardRateCoefB[r]) * exp(-backwardRateCoefC[r] / backwardReactionTemperature);
            
            Drf = (forwardRateCoefB[r] + forwardRateCoefC[r] / forwardReactionTemperature) / (forwardReactionTemperature * density * mixGasCv);
            Drb = (backwardRateCoefB[r] + backwardRateCoefC[r] / backwardReactionTemperature) / (backwardReactionTemperature * density * mixGasCv);

            //! Initialize.
            forwardRates[r] = 1.0;
            backwardRates[r] = 1.0;
            collisionRates[r] = 0.0;
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                densityPerMass = densityMassFraction[s] * speciesMassReciprocal[s];
                for (int i = 0; i < forwardStoichNumber[r][s]; ++ i)
                {
                    forwardRates[r] *= densityPerMass;
                }
                for (int i = 0; i < backwardStoichNumber[r][s]; ++ i)
                {
                    backwardRates[r] *= densityPerMass;
                }
                collisionRates[r] += densityPerMass * thirdBodyStoichNumber[r][s];
            }

            if (fabs(collisionRates[r]) < SMALL)     //! RCr.
            {
                collisionRates[r] = 1.0;
            }
            forwardRates[r] *= kf;    //! RAr.
            backwardRates[r] *= kb;   //! RBr.

            //! Compute the partial derivatives to density.
            alpha = 0.5 * squareVelocity - energy[numberOfSpecies - 1];
            theta[0] = forwardRates[r] * (forwardStoichNumber[r][numberOfSpecies - 1] / densityMassFraction[numberOfSpecies - 1] + Drf * alpha);
            theta[1] = backwardRates[r] * (backwardStoichNumber[r][numberOfSpecies - 1] / densityMassFraction[numberOfSpecies - 1] + Drb * alpha);
            theta[2] = thirdBodyStoichNumber[r][numberOfSpecies - 1] * speciesMassReciprocal[numberOfSpecies - 1];

            //! Compute the partial derivatives to density.
            partialValues[r][IR] = (theta[0] - theta[1]) * collisionRates[r] + (forwardRates[r] - backwardRates[r]) * theta[2];

            //! Compute the partial derivatives to rho*U, rho*V, rho*W and rho*E.
            phir = (Drf * forwardRates[r] - Drb * backwardRates[r]) * collisionRates[r];
            partialValues[r][IU] = - phir * uVelocity;
            partialValues[r][IV] = - phir * vVelocity;
            partialValues[r][IW] = - phir * wVelocity;
            partialValues[r][IP] = phir;

            //! Compute the partial derivatives to each species.
            for (int s = 0; s < numberOfSpecies - 1; ++ s)
            {
                gamaj = energy[numberOfSpecies - 1] - energy[s];
                theta[0] = forwardRates[r] * (forwardStoichNumber[r][s] / densityMassFraction[s] - forwardStoichNumber[r][numberOfSpecies - 1] / densityMassFraction[numberOfSpecies - 1] + Drf * gamaj);
                theta[1] = backwardRates[r] * (backwardStoichNumber[r][s] / densityMassFraction[s] - backwardStoichNumber[r][numberOfSpecies - 1] / densityMassFraction[numberOfSpecies - 1] + Drb * gamaj);
                theta[2] = thirdBodyStoichNumber[r][s] * speciesMassReciprocal[s] - thirdBodyStoichNumber[r][numberOfSpecies - 1] * speciesMassReciprocal[numberOfSpecies - 1];
                //! Compute the partial derivatives to the species j.
                partialValues[r][nm + s] = (theta[0] - theta[1]) * collisionRates[r] + (forwardRates[r] - backwardRates[r]) * theta[2];
            }
            //! Compute the partial derivative to vibration and electron energy.
            for (int i = 1; i < ntmodel; ++ i)
            {
                partialValues[r][nl + i] = - phir;
            }
        }
        //! Compute the chemical source term.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            chemicalSourceTerms[s][iCell] = 0.0;
            for (int r = 0; r < numberOfReaction; ++ r)
            {
                chemicalSourceTerms[s][iCell] += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * (forwardRates[r] - backwardRates[r]) * collisionRates[r];
            }
            chemicalSourceTerms[s][iCell] *= referenceSpeciesEquation * cellVolume[iCell] * speciesMass[s];
            //chemicalSourceTerms[s][iCell] *= speciesMass[s];
        }

        //! Construct elements of the Jacobian matrix of the chemical source term.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            //! Compute partial derivatives to species mass fraction.
            for (int j = 0; j < nEquation; ++ j)
            {
                sourceDerivatives[s][j][iCell] = 0.0;
                for (int r = 0; r < numberOfReaction; ++ r)
                {
                    sourceDerivatives[s][j][iCell] += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * partialValues[r][j];
                }
                sourceDerivatives[s][j][iCell] *= reference1 * cellVolume[iCell] * speciesMass[s];
                }
            //! Translate to the non-dimensional values.
            sourceDerivatives[s][IU][iCell] *= referenceVelocity;
            sourceDerivatives[s][IV][iCell] *= referenceVelocity;
            sourceDerivatives[s][IW][iCell] *= referenceVelocity;
            sourceDerivatives[s][IP][iCell] *= referenceVelocity * referenceVelocity;
            for (int i = 1; i < ntmodel; ++ i)
            {
                sourceDerivatives[s][nl + i][iCell] *= referenceVelocity * referenceVelocity;
            }
        }
    }

    delete []energy;
    delete []forwardRates;
    delete []backwardRates;
    delete []collisionRates;
    delete []massFraction;
    delete []densityMassFraction;
    DelPointer2(partialValues);
}

void MixedGas::ComputeChemicalSourceDerivatives2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{
    //! The function is used for chemical reaction including electron.

    using namespace IDX;     //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();

    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();

    //! The dimensional value of control temperature of chemical reaction.
    RDouble forwardReactionTemperature, backwardReactionTemperature;
    RDouble densityPerMass, density, squareVelocity, mixGasCv = 0.0;

    //! For source term in species transport equation.
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble reference1 = referenceDensity * referenceSpeciesEquation;

    RDouble *forwardRateCoefA = reactionRate->GetForwardReactionRateCoefficientA();
    RDouble *forwardRateCoefB = reactionRate->GetForwardReactionRateCoefficientB();
    RDouble *forwardRateCoefC = reactionRate->GetForwardReactionRateCoefficientC();

    RDouble *backwardRateCoefA = reactionRate->GetBackwardReactionRateCoefficientA();
    RDouble *backwardRateCoefB = reactionRate->GetBackwardReactionRateCoefficientB();
    RDouble *backwardRateCoefC = reactionRate->GetBackwardReactionRateCoefficientC();

    int **forwardStoichNumber = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
    int **backwardStoichNumber = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
    RDouble **thirdBodyStoichNumber = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();

    //! Record the coefficients in the chemical reactions.
    RDouble *forwardRates = new RDouble [numberOfReaction];
    RDouble *backwardRates = new RDouble [numberOfReaction];
    RDouble *collisionRates = new RDouble [numberOfReaction];
    RDouble *energy = new RDouble [numberOfSpecies];
    RDouble *massFraction = new RDouble [numberOfSpecies];
    RDouble *densityMassFraction = new RDouble [numberOfSpecies];
    int nEquation = this->nm + this->numberOfSpecies + this->ntmodel - 1;
    RDouble **partialValues = NewPointer2 <RDouble> (numberOfReaction, nEquation);

    //! Compute the elements of Jacobian matrix.
    RDouble coefA = 0.0, Drf = 0.0, Drb = 0.0, phir = 0.0, alpha = 0.0, gamaj = 0.0;
    RDouble uVelocity, vVelocity, wVelocity, theta[3] = {0.0};
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();
    for (int iCell = 0; iCell < nCellNumber; ++ iCell)
    {
        //! The control temperature of front reaction.
        forwardReactionTemperature = MAX(temperatures[ITT][iCell] * referenceTemperature, SMALL);
        //! The control temperature of back reaction. it is temperature of translation and rotation.
        backwardReactionTemperature = forwardReactionTemperature;

        density = primitiveVariables[IR ][iCell] * referenceDensity;
        uVelocity = primitiveVariables[IU][iCell] * referenceVelocity;
        vVelocity = primitiveVariables[IV][iCell] * referenceVelocity;
        wVelocity = primitiveVariables[IW][iCell] * referenceVelocity;
        squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            massFraction[s] = primitiveVariables[nm + s][iCell];
            densityMassFraction[s] = MAX(SMALL, density * massFraction[s]);
        }

        if (ntmodel == 1)     //! One-Temperature Model.
        {
            //! Compute specific heat at constant volume of each species.
            ComputeSpeciesConstantVolumeSpecificHeatDimensional(temperatures[ITT][iCell], energy);
            ComputeMixtureByMassFraction(massFraction, energy, mixGasCv);
            for (int s = 0; s < numberOfSpecies; ++ s)     //! Compute energy of each species.
            {
                energy[s] *= forwardReactionTemperature;
            }
        }
        else       //! Multi-Temperature Model.
        {
            //! Compute specific heat at constant volume of each species.
            ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(energy);
            ComputeMixtureByMassFraction(massFraction, energy, mixGasCv);
            for (int s = 0; s < numberOfSpecies; ++ s) // Compute energy of each species.
            {
                energy[s] = energy[s] * forwardReactionTemperature + thermodynamicProperties[s].hs0;
            }
        }

        for (int r = 0; r < numberOfReaction; ++ r)
        {
            //! Compute equilibrium constant of the front reaction.
            RDouble kf = forwardRateCoefA[r] * pow(forwardReactionTemperature, forwardRateCoefB[r]) * exp(-forwardRateCoefC[r] / forwardReactionTemperature);
            //! Compute equilibrium constant of the reverse reaction.
            RDouble kb = backwardRateCoefA[r] * pow(backwardReactionTemperature, backwardRateCoefB[r]) * exp(-backwardRateCoefC[r] / backwardReactionTemperature);

            Drf = (forwardRateCoefB[r] + forwardRateCoefC[r] / forwardReactionTemperature) / (forwardReactionTemperature * density * mixGasCv);
            Drb = (backwardRateCoefB[r] + backwardRateCoefC[r] / backwardReactionTemperature) / (backwardReactionTemperature * density * mixGasCv);

            //! Initialize.
            forwardRates[r] = 1.0;
            backwardRates[r] = 1.0;
            collisionRates[r] = 0.0;
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                densityPerMass = densityMassFraction[s] * speciesMassReciprocal[s];
                for (int i = 0; i < forwardStoichNumber[r][s]; ++ i)
                {
                    forwardRates[r] *= densityPerMass;
                }
                for (int i = 0; i < backwardStoichNumber[r][s]; ++ i)
                {
                    backwardRates[r] *= densityPerMass;
                }
                collisionRates[r] += densityPerMass * thirdBodyStoichNumber[r][s];
            }

            if (fabs(collisionRates[r]) < SMALL)    //! RCr.
            {
                collisionRates[r] = 1.0;
            }
            forwardRates[r] *= kf;    //! RAr.
            backwardRates[r] *= kb;   //! RBr.

            //! Compute the partial derivatives to density.
            alpha = 0.5 * squareVelocity - energy[numberOfSpecies - 2];
            theta[0] = forwardRates[r] * (forwardStoichNumber[r][numberOfSpecies - 2] / densityMassFraction[numberOfSpecies - 2] + Drf * alpha);
            theta[1] = backwardRates[r] * (backwardStoichNumber[r][numberOfSpecies - 2] / densityMassFraction[numberOfSpecies - 2] + Drb * alpha);
            theta[2] = thirdBodyStoichNumber[r][numberOfSpecies - 1] * speciesMassReciprocal[numberOfSpecies - 1];
            theta[2]+= thirdBodyStoichNumber[r][numberOfSpecies - 2] * speciesMassReciprocal[numberOfSpecies - 2];

            //! Compute the partial derivatives to density.
            partialValues[r][IR] = (theta[0] - theta[1]) * collisionRates[r] + (forwardRates[r] - backwardRates[r]) * theta[2];
            //! Compute the partial derivatives to rho*U, rho*V, rho*W and rho*E.
            phir = (Drf * forwardRates[r] - Drb * backwardRates[r]) * collisionRates[r];
            partialValues[r][IU] = - phir * uVelocity;
            partialValues[r][IV] = - phir * vVelocity;
            partialValues[r][IW] = - phir * wVelocity;
            partialValues[r][IP] = phir;

            //! Compute the partial derivatives to each species.
            for (int s = 0; s < numberOfSpecies - 2; ++ s)
            {
                coefA = ionType[s] * speciesMassReciprocal[s] / speciesMassReciprocal[numberOfSpecies - 1];
                gamaj = energy[numberOfSpecies - 2] - energy[s] + coefA * (energy[numberOfSpecies - 2] - energy[numberOfSpecies - 1]);

                theta[0]  = forwardStoichNumber[r][s] / densityMassFraction[s] - forwardStoichNumber[r][numberOfSpecies - 2] / densityMassFraction[numberOfSpecies - 2];
                theta[0] += coefA * (forwardStoichNumber[r][numberOfSpecies - 1] / densityMassFraction[numberOfSpecies - 1] - forwardStoichNumber[r][numberOfSpecies - 2] / densityMassFraction[numberOfSpecies - 2]);
                theta[0] += Drf * gamaj;
                theta[0]  = forwardRates[r] * theta[0];

                theta[1]  = backwardStoichNumber[r][s] / densityMassFraction[s] - backwardStoichNumber[r][numberOfSpecies - 2] / densityMassFraction[numberOfSpecies - 2];
                theta[1] += coefA * (backwardStoichNumber[r][numberOfSpecies - 1] / densityMassFraction[numberOfSpecies - 1] - backwardStoichNumber[r][numberOfSpecies - 2] / densityMassFraction[numberOfSpecies - 2]);
                theta[1] += Drb * gamaj;
                theta[1]  = backwardRates[r] * theta[1];

                theta[2]  = thirdBodyStoichNumber[r][s] * speciesMassReciprocal[s] - thirdBodyStoichNumber[r][numberOfSpecies - 2] * speciesMassReciprocal[numberOfSpecies - 2];
                theta[2] += coefA * (thirdBodyStoichNumber[r][numberOfSpecies - 1] * speciesMassReciprocal[numberOfSpecies - 1] - thirdBodyStoichNumber[r][numberOfSpecies - 2] * speciesMassReciprocal[numberOfSpecies - 2]);

                partialValues[r][nm + s] = (theta[0] - theta[1]) * collisionRates[r] + (forwardRates[r] - backwardRates[r]) * theta[2];
            }
            //! Compute the partial derivative to vibration and electron energy.
            for (int i = 1; i < ntmodel; ++ i)
            {
                partialValues[r][nl + i] = - phir;
            }
        }
        //! Compute the chemical source term.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            chemicalSourceTerms[s][iCell] = 0.0;
            for (int r = 0; r < numberOfReaction; ++ r)
            {
                chemicalSourceTerms[s][iCell] += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * (forwardRates[r] - backwardRates[r]) * collisionRates[r];
            }
            chemicalSourceTerms[s][iCell] *= referenceSpeciesEquation * cellVolume[iCell] * speciesMass[s];
            //chemicalSourceTerms[s][iCell] *= speciesMass[s];
        }

        //! Construct elements of the Jacobian matrix of the chemical source term.
        for (int s = 0; s < numberOfSpecies - 1; ++ s)
        {
            //! Compute partial derivatives to species mass fraction.
            for (int j = 0; j < nEquation; ++ j)
            {
                sourceDerivatives[s][j][iCell] = 0.0;
                for (int r = 0; r < numberOfReaction; ++ r)
                {
                    sourceDerivatives[s][j][iCell] += (backwardStoichNumber[r][s] - forwardStoichNumber[r][s]) * partialValues[r][j];
                }
                sourceDerivatives[s][j][iCell] *= reference1 * cellVolume[iCell] * speciesMass[s];
            }
            //! Translate to the non-dimensional values.
            sourceDerivatives[s][IU][iCell] *= referenceVelocity;
            sourceDerivatives[s][IV][iCell] *= referenceVelocity;
            sourceDerivatives[s][IW][iCell] *= referenceVelocity;
            sourceDerivatives[s][IP][iCell] *= referenceVelocity * referenceVelocity;
            for (int i = 1; i < ntmodel; ++ i)
            {
                sourceDerivatives[s][nl + i][iCell] *= referenceVelocity * referenceVelocity;
            }
        }
    }

    delete []energy;
    delete []forwardRates;
    delete []backwardRates;
    delete []collisionRates;
    delete []massFraction;
    delete []densityMassFraction;
    DelPointer2(partialValues);
}

void MixedGas::ComputeEnergySourceDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                              RDouble ***chemicalSourceDerivatives, RDouble ***energySourceDerivatives, int nCellNumber)
{
    //! The function is used for two-temperature model.

    using namespace IDX;     //! The following variables is dimensional values.
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermoProperty = molecularProperty->GetThermodynamicProperties();
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();

    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble referenceOmega = reynoldsReferenceLengthDimensional / referenceVelocity;
    RDouble reference1 = reynoldsReferenceLengthDimensional / pow(referenceVelocity, 3);

    RDouble *Omegas = new RDouble[numberOfSpecies];
    int nEquation = nm + numberOfSpecies + ntmodel - 1;
    RDouble **omegaDerivatives = NewPointer2 <RDouble> (numberOfSpecies, nEquation);
    RDouble **energyDerivatives = NewPointer2 <RDouble> (ntmodel - 1, nEquation);
    RDouble *massFraction = new RDouble[numberOfSpecies];
    RDouble *speciesCvtr = new RDouble[numberOfSpecies];
    RDouble *speciesCvv = new RDouble[numberOfSpecies];
    RDouble *speciesCve = new RDouble[numberOfSpecies];
    RDouble *speciesCvvStar1 = new RDouble[numberOfSpecies];
    RDouble *speciesCvvStar2 = new RDouble[numberOfSpecies];
    RDouble *speciesTransRotationEnergy = new RDouble[numberOfSpecies];
    RDouble *speciesVibrationEnergy = new RDouble[numberOfSpecies];
    RDouble *speciesElectronEnergy = new RDouble[numberOfSpecies];
    RDouble *speciesVibrationEnergyStar1 = new RDouble[numberOfSpecies];
    RDouble *speciesVibrationEnergyStar2 = new RDouble[numberOfSpecies];
    RDouble *speciesVibrationTime = new RDouble[numberOfSpecies];
    RDouble *speciesElectronTime = new RDouble[numberOfSpecies];

    //! Compute the partial derivatives.
    RDouble mixGasCvtr, mixGasCvv, mixGasCve, transRotationTemperature, vibrationTemperature, electronTemperature;
    RDouble coefA = 0.0, energy = 0.0, energyj = 0.0, electronPressure = 0.0;
    RDouble gasConstant = this->GetUniversalGasConstant();
    for (int iCell = 0; iCell < nCellNumber; ++ iCell)
    {
        //! Obtain the primitive variables.
        RDouble density = primitiveVariables[IR][iCell] * referenceDensity;
        RDouble uVelocity = primitiveVariables[IU][iCell] * referenceVelocity;
        RDouble vVelocity = primitiveVariables[IV][iCell] * referenceVelocity;
        RDouble wVelocity = primitiveVariables[IW][iCell] * referenceVelocity;
        RDouble pressure = primitiveVariables[IP][iCell] * referenceDensity * referenceVelocity * referenceVelocity;
        RDouble squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            //! Obtain the chemical source terms.
            Omegas[s] = chemicalSourceTerms[s][iCell] / (referenceSpeciesEquation * cellVolume[iCell]);
            //! Obtain the partial derivatives of chemical source to the vibration and electron energy.
            for (int i = 0; i < nEquation; ++ i)
            {
                omegaDerivatives[s][i] = chemicalSourceDerivatives[s][i][iCell] / (referenceOmega * cellVolume[iCell]);
            }
            omegaDerivatives[s][IU] /= referenceVelocity;
            omegaDerivatives[s][IV] /= referenceVelocity;
            omegaDerivatives[s][IW] /= referenceVelocity;
            omegaDerivatives[s][IP] /= referenceVelocity * referenceVelocity;
            for (int i = 1; i < ntmodel; ++ i)
            {
                omegaDerivatives[s][nl + i] /= referenceVelocity * referenceVelocity;
            }
        }

        //! Obtain the mass fraction of species.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            massFraction[s] = primitiveVariables[nm + s][iCell];
        }
        //! Obtain the temperature of two-temperature model.
        transRotationTemperature = temperatures[ITT][iCell];
        vibrationTemperature = temperatures[ITV][iCell];
        electronTemperature = vibrationTemperature;

        //! Compute specific heat at constant volume.
        this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCvtr);
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(vibrationTemperature, speciesCvv);
        this->ComputeDimensionalElectronSpecificHeatAtConstantVolume(electronTemperature, speciesCve);
        this->ComputeMixtureByMassFraction(massFraction, speciesCvtr, mixGasCvtr);
        this->ComputeMixtureByMassFraction(massFraction, speciesCvv, mixGasCvv);
        this->ComputeMixtureByMassFraction(massFraction, speciesCve, mixGasCve);
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(transRotationTemperature, speciesCvvStar1);

        //! Compute the energies of species.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            speciesTransRotationEnergy[s] = speciesCvtr[s] * transRotationTemperature * referenceTemperature + thermoProperty[s].hs0;
        }
        this->ComputeDimensionalVibrationEnergy(vibrationTemperature, speciesVibrationEnergy);
        this->ComputeDimensionalElectronEnergy(electronTemperature, speciesElectronEnergy);
        //! Compute the vibration energy by the translation-rotation temperature.
        this->ComputeDimensionalVibrationEnergy(transRotationTemperature, speciesVibrationEnergyStar1);

        //! Compute the relaxation time.of translation-vibration mode.
        this->GetTranslationAndVibrationRelaxationTime(thermoProperty, massFraction, speciesMass, transRotationTemperature * referenceTemperature,
                                                       pressure, density, speciesVibrationTime);
        RDouble part1 = 0.0, part2 = 0.0, part3 = 0.0;
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            if (thermoProperty[s].nType < 2)
            {
                continue;
            }
            part1 += massFraction[s] * speciesCvvStar1[s] / speciesVibrationTime[s];
            part2 += Omegas[s] * (speciesCvv[s] + speciesCve[s]);
            part3 += massFraction[s] * speciesCvv[s] / speciesVibrationTime[s];
        }
        part1 /= mixGasCvtr;
        part2 /= density * (mixGasCvv + mixGasCve) + SMALL;
        part3 /= mixGasCvv + mixGasCve + SMALL;

        //! Compute the partial derivatives of c-ve term.
        for (int i = 0; i < nEquation; ++ i)
        {
            energyDerivatives[0][i] = 0.0;
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                if (thermoProperty[s].nType < 2)
                {
                    continue;
                }
                energyDerivatives[0][i] += omegaDerivatives[s][i] * (speciesVibrationEnergy[s] + speciesElectronEnergy[s]);
            }
        }

        //! Compute the partial derivatives of t-v term.
        if (this->nElectronIndex < 0)      //! No electron.
        {
            //! Partial derivative to density.
            energy = 0.5 * squareVelocity - speciesTransRotationEnergy[numberOfSpecies - 1];
            RDouble temp = (speciesVibrationEnergyStar1[numberOfSpecies - 1] - speciesVibrationEnergy[numberOfSpecies - 1]) / speciesVibrationTime[numberOfSpecies - 1];
            energyDerivatives[0][IR] += energy * part1 + temp;
            //! Partial derivative to density * velocity.
            energyDerivatives[0][IU] -= uVelocity * part1;
            energyDerivatives[0][IV] -= vVelocity * part1;
            energyDerivatives[0][IW] -= wVelocity * part1;
            energyDerivatives[0][IP] += part1;
            //! Partial derivative to density * species.
            for (int s = 0; s < numberOfSpecies - 1; ++ s)
            {
                energyj = speciesTransRotationEnergy[numberOfSpecies - 1] - speciesTransRotationEnergy[s];
                energyDerivatives[0][nm + s] += energyj * part1;
                energyDerivatives[0][nm + s] += (speciesVibrationEnergyStar1[s] - speciesVibrationEnergy[s]) / speciesVibrationTime[s] - temp;
            }
            //! Partial derivative to vibration energy.
            energyDerivatives[0][nl + 1] += part2 - part1 - part3;
        }
        else        //! The species including the electron.
        {
            //! Compute the vibration energy by the electron temperature.
            this->ComputeDimensionalVibrationEnergy(electronTemperature, speciesVibrationEnergyStar2);

            //! Compute the relaxation time of electron-vibration mode.
            electronPressure = density * massFraction[numberOfSpecies - 1] * gasConstant * speciesMassReciprocal[numberOfSpecies - 1] * electronTemperature * referenceTemperature;
            this->GetElectronAndVibrationRelaxationTime(thermoProperty, electronPressure, electronTemperature * referenceTemperature, speciesElectronTime);
        }

        //! To obtain the non-dimensional values.
        for (int j = 0; j < nEquation; ++ j)
        {
            energySourceDerivatives[0][j][iCell] = reference1 * energyDerivatives[0][j] * cellVolume[iCell];
        }
        energySourceDerivatives[0][IU][iCell] *= referenceVelocity;
        energySourceDerivatives[0][IV][iCell] *= referenceVelocity;
        energySourceDerivatives[0][IW][iCell] *= referenceVelocity;
        energySourceDerivatives[0][IP][iCell] *= referenceVelocity * referenceVelocity;
        energySourceDerivatives[0][nl + 1][iCell] *= referenceVelocity * referenceVelocity;
    }

    delete []Omegas;
    delete []massFraction;
    delete []speciesCvtr;
    delete []speciesCvv;
    delete []speciesCve;
    delete []speciesCvvStar1;
    delete []speciesCvvStar2;
    delete []speciesTransRotationEnergy;
    delete []speciesVibrationEnergy;
    delete []speciesElectronEnergy;
    delete []speciesVibrationEnergyStar1;
    delete []speciesVibrationEnergyStar2;
    delete []speciesVibrationTime;
    delete []speciesElectronTime;
    DelPointer2(omegaDerivatives);
    DelPointer2(energyDerivatives);
}

void MixedGas::ComputeEnergySourceDiagonalDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                                      RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber)
{
    //! The function is used for two-temperature model.

    using namespace IDX;    //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble referenceEnergyEquation = referenceVelocity * reynoldsReferenceLengthDimensional;
    RDouble reference1 = reynoldsReferenceLengthDimensional / referenceVelocity;
    RDouble reference2 = referenceSpeciesEquation / pow(referenceVelocity, 2);

    RDouble *massFraction = new RDouble[numberOfSpecies];
    RDouble *Omegas = new RDouble[numberOfSpecies];
    RDouble *sigma = new RDouble[numberOfSpecies];
    RDouble *speciesCvv = new RDouble[numberOfSpecies];
    RDouble *speciesCve = new RDouble[numberOfSpecies];
    RDouble *speciesCvvStar1 = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergy = new RDouble[numberOfSpecies];
    RDouble *electronEnergy = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergyStar1 = new RDouble[numberOfSpecies];
    RDouble *transVibrationTime = new RDouble[numberOfSpecies];
    RDouble **OmegaDerivative = NewPointer2 <RDouble> (numberOfSpecies, ntmodel - 1);

    //! Compute the partial derivatives in the main diagonal position.
    RDouble mixGasCvtr, mixGasCvv, mixGasCve, mixGasCvve;
    RDouble temperature[3], Ttr = 0.0, Tve = 0.0, vibrationOmega = 0.0;
    RDouble gasConstant = rjmk, density = 0.0, pressure = 0.0;
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermoProperty = molecularProperty->GetThermodynamicProperties();

    for (int iCell = 0; iCell < nCellNumber; ++ iCell)
    {
        density = primitiveVariables[IR][iCell] * referenceDensity;
        pressure = primitiveVariables[IP][iCell] * referenceDensity * referenceVelocity * referenceVelocity;

        for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
        {
            //! Obtain the partial derivatives of chemical source to the vibration and electron energy.
            Omegas[s] = chemicalSourceTerms[s][iCell] / (referenceSpeciesEquation * cellVolume[iCell]);
            for (int it = 1; it < ntmodel; ++ it)
            {
                OmegaDerivative[s][it - 1] = chemicalSourceDerivatives[s][nl + it][iCell] / (referenceEnergyEquation * cellVolume[iCell]);
            }
        }
        //! Obtain the mass fraction of species.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            massFraction[s] = primitiveVariables[nm + s][iCell];
        }
        //! Obtain the temperature.
        for (int it = 0; it < ntmodel; ++ it)
        {
            temperature[it] = temperatures[it][iCell];
        }
        temperature[ITE] = temperature[ITV];
        Ttr = temperature[ITT] * referenceTemperature;
        Tve = temperature[ITV] * referenceTemperature;

        //! Compute the specific heat at constant volume.
        mixGasCvtr = this->GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(massFraction);
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITV], speciesCvv);
        this->ComputeDimensionalElectronSpecificHeatAtConstantVolume(temperature[ITE], speciesCve);
        this->ComputeMixtureByMassFraction(massFraction, speciesCvv, mixGasCvv);
        this->ComputeMixtureByMassFraction(massFraction, speciesCve, mixGasCve);
        mixGasCvve = mixGasCvv + mixGasCve;

        //! Compute vibration energy by translation-rotation temperature.
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITT], speciesCvvStar1);
        //! Compute vibration energy and electron energy.
        this->ComputeDimensionalVibrationEnergy(temperature[ITV], vibrationEnergy);
        this->ComputeDimensionalElectronEnergy(temperature[ITE], electronEnergy);
        this->ComputeDimensionalVibrationEnergy(temperature[ITT], vibrationEnergyStar1);

        //! Compute the relaxation time.of translation-vibration mode.
        this->GetTranslationAndVibrationRelaxationTime(thermoProperty, massFraction, speciesMass, Ttr, pressure, density, transVibrationTime);

        RDouble part1 = 0.0, part2 = 0.0, part3 = 0.0, part4 = 0.0;
        vibrationOmega = 0.0;
        for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
        {
            if (thermoProperty[s].nType < 2)
            {
                continue;
            }

            part1 += massFraction[s] * speciesCvvStar1[s] / transVibrationTime[s];
            part2 += massFraction[s] * speciesCvv[s] / transVibrationTime[s];
            part3 += Omegas[s] * (speciesCvv[s] + speciesCve[s]);
            part4 += OmegaDerivative[s][0] * (vibrationEnergy[s] + electronEnergy[s]);

            //! Compute the source term of vibration energy equation.
            vibrationOmega += density * massFraction[s] * (vibrationEnergyStar1[s] - vibrationEnergy[s]) / transVibrationTime[s];
            vibrationOmega += Omegas[s] * (vibrationEnergy[s] + electronEnergy[s]);
        }
        part1 = part1 / mixGasCvtr;
        part2 = part2 / mixGasCvve;
        part3 = part3 / (density * mixGasCvve);
        //! Compute derivatives.
        energySourceDerivatives[0][iCell] = - part1 - part2 + part3 + part4;

        if (this->nElectronIndex >= 0) //species including electron.
        {
            //! Compute the number density of electron.
            RDouble electronDensity = density * massFraction[nElectronIndex];
            RDouble Ne = AVOGADRO_CONSTANT * electronDensity * speciesMassReciprocal[nElectronIndex];
            this->GetEffectiveCollisionCrossSectionalArea(thermoProperty, Tve, Ne, sigma);

            RDouble Fe = 0.0, Ns = 0.0, fes = 0.0, speciesDensity = 0.0;
            //! The heat motion velocity of electron.
            RDouble heatVelocity = sqrt(8.0 * gasConstant * Tve * speciesMassReciprocal[nElectronIndex] / PI);
            for (int s = 0; s <= nLastSpeciesIndex; ++ s)
            {
                speciesDensity = density * massFraction[s];
                Ns = AVOGADRO_CONSTANT * speciesDensity * speciesMassReciprocal[s];
                fes = Ns * sigma[s] * heatVelocity;
                Fe += fes * speciesMassReciprocal[s];
            }
            //! Compute the source term of vibration energy equation.
            vibrationOmega += 3.0 * electronDensity * gasConstant * (Ttr - Tve);

            RDouble mixGasEv = 0.0, mixGasEe = 0.0, mixGasEve;
            this->ComputeMixtureByMassFraction(massFraction, vibrationEnergy, mixGasEv);
            this->ComputeMixtureByMassFraction(massFraction, electronEnergy , mixGasEe);
            mixGasEve = mixGasEv + mixGasEe;
            energySourceDerivatives[0][iCell] -= 3.0 * massFraction[nElectronIndex] * gasConstant * Fe * (1.0 / mixGasCvtr + 1.0 / mixGasCvve + (Tve - Ttr) / (2.0 * mixGasEve));
        }
        energySourceDerivatives[0][iCell] *= reference1 * cellVolume[iCell];
        //! The source term of vibration energy equation.
        chemicalSourceTerms[numberOfSpecies][iCell] = reference2 * vibrationOmega * cellVolume[iCell];
    }

    delete []Omegas;
    delete []sigma;
    delete []massFraction;
    delete []speciesCvv;
    delete []speciesCve;
    delete []speciesCvvStar1;
    delete []vibrationEnergy;
    delete []electronEnergy;
    delete []vibrationEnergyStar1;
    delete []transVibrationTime;
    DelPointer2(OmegaDerivative);
}

void MixedGas::ComputeEnergySourceDiagonalDerivatives2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                                      RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber)
{
    //! The function is used for three-temperature model.

    using namespace IDX;     //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble referenceEnergyEquation = referenceVelocity * reynoldsReferenceLengthDimensional;
    RDouble reference1 = reynoldsReferenceLengthDimensional / referenceVelocity;
    RDouble reference2 = referenceSpeciesEquation / pow(referenceVelocity, 2);

    RDouble *massFraction = new RDouble[numberOfSpecies];
    RDouble *Omegas = new RDouble[numberOfSpecies];
    RDouble *sigma = new RDouble[numberOfSpecies];
    RDouble *speciesCvv = new RDouble[numberOfSpecies];
    RDouble *speciesCve = new RDouble[numberOfSpecies];
    RDouble *speciesCvvStar1 = new RDouble[numberOfSpecies];
    RDouble *speciesCvvStar2 = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergy = new RDouble[numberOfSpecies];
    RDouble *electronEnergy = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergyStar1 = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergyStar2 = new RDouble[numberOfSpecies];
    RDouble *transVibrationTime = new RDouble[numberOfSpecies];
    RDouble *transElectronTime = new RDouble[numberOfSpecies];
    RDouble **OmegaDerivative = NewPointer2 <RDouble> (numberOfSpecies, ntmodel - 1);

    //! Compute the partial derivatives in the main diagonal position.
    RDouble gasConstant = rjmk, mixGasCvtr, mixGasCvv, mixGasCve;
    RDouble temperature[3], Ttr = 0.0, Te = 0.0, electronPressure = 0.0, electronDensity = 0.0;
    RDouble vibrationOmega = 0.0, electronOmega = 0.0;
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermoProperty = molecularProperty->GetThermodynamicProperties();

    for (int iCell = 0; iCell < nCellNumber; ++ iCell)
    {
        RDouble density = primitiveVariables[IR][iCell] * referenceDensity;
        RDouble pressure = primitiveVariables[IP][iCell] * referenceDensity * referenceVelocity * referenceVelocity;

        for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
        {
            //! Obtain the partial derivatives of chemical source to the vibration and electron energy.
            Omegas[s] = chemicalSourceTerms[s][iCell] / (referenceSpeciesEquation * cellVolume[iCell]);
            for (int it = 1; it < ntmodel; ++ it)
            {
                OmegaDerivative[s][it - 1] = chemicalSourceDerivatives[s][nl + it][iCell] / (referenceEnergyEquation * cellVolume[iCell]);
            }
        }
        //! Obtain the mass fraction of species.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            massFraction[s] = primitiveVariables[nm + s][iCell];
        }
        //! Obtain the temperature.
        for (int it = 0; it < ntmodel; ++ it)
        {
            temperature[it] = temperatures[it][iCell];
        }
        Ttr = temperature[ITT] * referenceTemperature;
        Te  = temperature[ITE] * referenceTemperature;

        //! Compute the specific heat at constant volume.
        mixGasCvtr = this->GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(massFraction);
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITV], speciesCvv);
        this->ComputeDimensionalElectronSpecificHeatAtConstantVolume(temperature[ITE], speciesCve);
        this->ComputeMixtureByMassFraction(massFraction, speciesCvv, mixGasCvv);
        this->ComputeMixtureByMassFraction(massFraction, speciesCve, mixGasCve);

        //! Compute vibration energy by translation-rotation temperature.
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITT], speciesCvvStar1);
        //! Compute vibration energy by electron temperature.
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITE], speciesCvvStar2);

        //! Compute vibration energy and electron energy.
        this->ComputeDimensionalVibrationEnergy(temperature[ITV], vibrationEnergy);
        this->ComputeDimensionalElectronEnergy(temperature[ITE], electronEnergy);
        this->ComputeDimensionalVibrationEnergy(temperature[ITT], vibrationEnergyStar1);
        this->ComputeDimensionalVibrationEnergy(temperature[ITE], vibrationEnergyStar2);

        //! Compute the relaxation time.of translation-vibration mode.
        this->GetTranslationAndVibrationRelaxationTime(thermoProperty, massFraction, speciesMass, Ttr, pressure, density, transVibrationTime);
        if (this->nElectronIndex >= 0)
        {
            electronDensity = density * massFraction[this->nElectronIndex];
            electronPressure = electronDensity * gasConstant * speciesMassReciprocal[this->nElectronIndex] * Te;
        }
        //! Compute the relaxation time.of electron-vibration mode.
        this->GetElectronAndVibrationRelaxationTime(thermoProperty, electronPressure, Te, transElectronTime);

        RDouble partv1 = 0.0, partv2 = 0.0, partv3 = 0.0, partv4 = 0.0;
        RDouble parte1 = 0.0, parte2 = 0.0, parte3 = 0.0;
        vibrationOmega = 0.0;
        electronOmega  = 0.0;
        for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
        {
            parte1 += Omegas[s] * speciesCve[s];
            parte2 += OmegaDerivative[s][1] * electronEnergy[s];
            //! Compute source term of electron energy equation.
            electronOmega += Omegas[s] * electronEnergy[s];

            if (thermoProperty[s].nType < 2)
            {
                continue;
            }
            partv1 += massFraction[s] * speciesCvvStar1[s] / transVibrationTime[s];
            partv2 += massFraction[s] * speciesCvv[s] / transVibrationTime[s];
            partv3 += Omegas[s] * speciesCvv[s];
            partv4 += OmegaDerivative[s][0] * vibrationEnergy[s];

            //! Compute source term of vibration energy equation.
            vibrationOmega += density * massFraction[s] * (vibrationEnergyStar1[s] - vibrationEnergy[s]) / transVibrationTime[s];
            vibrationOmega += Omegas[s] * vibrationEnergy[s];
        }
        if (nNitrogenIndex >= 0)
        {
            partv2 += massFraction[nNitrogenIndex] * speciesCvv[nNitrogenIndex] / transElectronTime[nNitrogenIndex];
            parte3 += massFraction[nNitrogenIndex] * speciesCvvStar2[nNitrogenIndex] / transElectronTime[nNitrogenIndex];
            //! Compute source term.
            RDouble temp = density * massFraction[nNitrogenIndex] * (vibrationEnergyStar2[nNitrogenIndex] - vibrationEnergy[nNitrogenIndex]) / transElectronTime[nNitrogenIndex];
            vibrationOmega += temp;
            electronOmega -= temp;
        }

        partv1 = partv1 / mixGasCvtr;
        partv2 = partv2 / mixGasCvv;
        partv3 = partv3 / (density * mixGasCvv);
        parte1 = parte1 / (density * mixGasCve);
        parte3 = parte3 / mixGasCve;

        //! Compute derivatives.
        energySourceDerivatives[0][iCell] = partv3 + partv4 - partv1 - partv2;
        energySourceDerivatives[1][iCell] = parte1 + parte2 - parte3;

        if (this->nElectronIndex >= 0)     //! Species including electron.
        {
            //! Compute the number density of electron.
            RDouble Ne = AVOGADRO_CONSTANT * electronDensity * speciesMassReciprocal[nElectronIndex];
            this->GetEffectiveCollisionCrossSectionalArea(thermoProperty, Te, Ne, sigma);

            RDouble Fe = 0.0, Ns = 0.0, fes = 0.0, speciesDensity = 0.0;
            RDouble heatVelocity = sqrt(8.0 * gasConstant * Te * speciesMassReciprocal[nElectronIndex] / PI);
            for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
            {
                speciesDensity = density * massFraction[s];
                Ns = AVOGADRO_CONSTANT * speciesDensity * speciesMassReciprocal[s];
                fes = Ns * sigma[s] * heatVelocity;
                Fe += fes * speciesMassReciprocal[s];
            }
            //! Compute source term of electron energy equation.
            electronOmega += 3.0 * electronDensity * gasConstant * (Ttr - Te) * Fe;

            RDouble mixGasEe = 0.0;
            this->ComputeMixtureByMassFraction(massFraction, electronEnergy , mixGasEe);
            energySourceDerivatives[1][iCell] -= 3.0 * massFraction[nElectronIndex] * gasConstant * Fe * (1.0 / mixGasCvtr + 1.0 / mixGasCve + (Te - Ttr) / (2.0 * mixGasEe));
        }
        energySourceDerivatives[0][iCell] *= reference1 * cellVolume[iCell];
        energySourceDerivatives[1][iCell] *= reference1 * cellVolume[iCell];

        //! The source terms of vibration energy equation and electron energy equation.
        chemicalSourceTerms[numberOfSpecies][iCell] = reference2 * vibrationOmega * cellVolume[iCell];
        chemicalSourceTerms[numberOfSpecies + 1][iCell] = reference2 * electronOmega * cellVolume[iCell];
    }

    delete []Omegas;
    delete []sigma;
    delete []massFraction;
    delete []speciesCvv;
    delete []speciesCve;
    delete []speciesCvvStar1;
    delete []speciesCvvStar2;
    delete []vibrationEnergy;
    delete []electronEnergy;
    delete []vibrationEnergyStar1;
    delete []vibrationEnergyStar2;
    delete []transVibrationTime;
    delete []transElectronTime;
    DelPointer2(OmegaDerivative);
}

void MixedGas::ComputeEnergySourceJacobianT2(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives)
{
    using namespace IDX;    //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble referenceEnergyEquation = referenceVelocity * reynoldsReferenceLengthDimensional;
    RDouble reference1 = reynoldsReferenceLengthDimensional / referenceVelocity;
    RDouble reference2 = referenceSpeciesEquation / pow(referenceVelocity, 2);

    RDouble *Omegas = new RDouble[numberOfSpecies];
    RDouble *sigma = new RDouble[numberOfSpecies];
    RDouble *speciesCvv = new RDouble[numberOfSpecies];
    RDouble *speciesCve = new RDouble[numberOfSpecies];
    RDouble *speciesCvvStar1 = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergy = new RDouble[numberOfSpecies];
    RDouble *electronEnergy = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergyStar1 = new RDouble[numberOfSpecies];
    RDouble *transVibrationTime = new RDouble[numberOfSpecies];
    RDouble **OmegaDerivative = NewPointer2 <RDouble> (numberOfSpecies, ntmodel - 1);

    //! Compute the partial derivatives in the main diagonal position.
    RDouble mixGasCvtr, mixGasCvv, mixGasCve, mixGasCvve;
    RDouble Ttr = 0.0, Tve = 0.0, vibrationOmega = 0.0;
    RDouble gasConstant = rjmk, density = 0.0, pressure = 0.0;
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermoProperty = molecularProperty->GetThermodynamicProperties();

    density = primitiveVariables[IR] * referenceDensity;
    pressure = primitiveVariables[IP] * referenceDensity * referenceVelocity * referenceVelocity;

    for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
    {
        //! Obtain the partial derivatives of chemical source to the vibration and electron energy.
        Omegas[s] = chemicalSourceTerms[s] / (referenceSpeciesEquation * cellVolume);
        for (int it = 1; it < ntmodel; ++ it)
        {
            OmegaDerivative[s][it - 1] = sourceDerivatives[s][numberOfSpecies + it - 1] / (referenceEnergyEquation * cellVolume);
        }
    }
    //! Obtain the mass fraction of species.
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        massFractions[s] = primitiveVariables[nm + s];
    }
    Ttr = temperatures[ITT] * referenceTemperature;
    Tve = temperatures[ITV] * referenceTemperature;

    //! Compute the specific heat at constant volume.
    mixGasCvtr = this->GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(massFractions);
    this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperatures[ITV], speciesCvv);
    this->ComputeDimensionalElectronSpecificHeatAtConstantVolume(temperatures[ITE], speciesCve);
    this->ComputeMixtureByMassFraction(massFractions, speciesCvv, mixGasCvv);
    this->ComputeMixtureByMassFraction(massFractions, speciesCve, mixGasCve);
    mixGasCvve = mixGasCvv + mixGasCve;

    //! Compute vibration energy by translation-rotation temperature.
    this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperatures[ITT], speciesCvvStar1);
    //! Compute vibration energy and electron energy.
    this->ComputeDimensionalVibrationEnergy(temperatures[ITV], vibrationEnergy);
    this->ComputeDimensionalElectronEnergy(temperatures[ITE], electronEnergy);
    this->ComputeDimensionalVibrationEnergy(temperatures[ITT], vibrationEnergyStar1);

    //! Compute the relaxation time.of translation-vibration mode.
    this->GetTranslationAndVibrationRelaxationTime(thermoProperty, massFractions, speciesMass, Ttr, pressure, density, transVibrationTime);

    RDouble part1 = 0.0, part2 = 0.0, part3 = 0.0, part4 = 0.0;
    vibrationOmega = 0.0;
    for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
    {
        if (thermoProperty[s].nType < 2)
        {
            continue;
        }

        part1 += massFractions[s] * speciesCvvStar1[s] / transVibrationTime[s];
        part2 += massFractions[s] * speciesCvv[s] / transVibrationTime[s];
        part3 += Omegas[s] * (speciesCvv[s] + speciesCve[s]);
        part4 += OmegaDerivative[s][0] * (vibrationEnergy[s] + electronEnergy[s]);

        //! Compute the source term of vibration energy equation.
        vibrationOmega += density * massFractions[s] * (vibrationEnergyStar1[s] - vibrationEnergy[s]) / transVibrationTime[s];
        vibrationOmega += Omegas[s] * (vibrationEnergy[s] + electronEnergy[s]);
    }
    part1 = part1 / mixGasCvtr;
    part2 = part2 / mixGasCvve;
    part3 = part3 / (density * mixGasCvve);
    //! Compute derivatives.
    sourceDerivatives[numberOfSpecies][numberOfSpecies] = - part1 - part2 + part3 + part4;

    if (this->nElectronIndex >= 0) //species including electron.
    {
        //! Compute the number density of electron.
        RDouble electronDensity = density * massFractions[nElectronIndex];
        RDouble Ne = AVOGADRO_CONSTANT * electronDensity * speciesMassReciprocal[nElectronIndex];
        this->GetEffectiveCollisionCrossSectionalArea(thermoProperty, Tve, Ne, sigma);

        RDouble Fe = 0.0, Ns = 0.0, fes = 0.0, speciesDensity = 0.0;
        //! The heat motion velocity of electron.
        RDouble heatVelocity = sqrt(8.0 * gasConstant * Tve * speciesMassReciprocal[nElectronIndex] / PI);
        for (int s = 0; s <= nLastSpeciesIndex; ++ s)
        {
            speciesDensity = density * massFractions[s];
            Ns = AVOGADRO_CONSTANT * speciesDensity * speciesMassReciprocal[s];
            fes = Ns * sigma[s] * heatVelocity;
            Fe += fes * speciesMassReciprocal[s];
        }
        //! Compute the source term of vibration energy equation.
        vibrationOmega += 3.0 * electronDensity * gasConstant * (Ttr - Tve);

        RDouble mixGasEv = 0.0, mixGasEe = 0.0, mixGasEve;
        this->ComputeMixtureByMassFraction(massFractions, vibrationEnergy, mixGasEv);
        this->ComputeMixtureByMassFraction(massFractions, electronEnergy , mixGasEe);
        mixGasEve = mixGasEv + mixGasEe;
        sourceDerivatives[numberOfSpecies][numberOfSpecies] -= 3.0 * massFractions[nElectronIndex] * gasConstant * Fe * (1.0 / mixGasCvtr + 1.0 / mixGasCvve + (Tve - Ttr) / (2.0 * mixGasEve));
    }
    sourceDerivatives[numberOfSpecies][numberOfSpecies] *= reference1 * cellVolume;
    //! The source term of vibration energy equation.
    chemicalSourceTerms[numberOfSpecies] = reference2 * vibrationOmega * cellVolume;

    delete []Omegas;
    delete []sigma;
    delete []speciesCvv;
    delete []speciesCve;
    delete []speciesCvvStar1;
    delete []vibrationEnergy;
    delete []electronEnergy;
    delete []vibrationEnergyStar1;
    delete []transVibrationTime;
    DelPointer2(OmegaDerivative);
}

void MixedGas::ComputeEnergySourceJacobianT3(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives)
{
    using namespace IDX;     //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble referenceEnergyEquation = referenceVelocity * reynoldsReferenceLengthDimensional;
    RDouble reference1 = reynoldsReferenceLengthDimensional / referenceVelocity;
    RDouble reference2 = referenceSpeciesEquation / pow(referenceVelocity, 2);

    RDouble *Omegas = new RDouble[numberOfSpecies];
    RDouble *sigma = new RDouble[numberOfSpecies];
    RDouble *speciesCvv = new RDouble[numberOfSpecies];
    RDouble *speciesCve = new RDouble[numberOfSpecies];
    RDouble *speciesCvvStar1 = new RDouble[numberOfSpecies];
    RDouble *speciesCvvStar2 = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergy = new RDouble[numberOfSpecies];
    RDouble *electronEnergy = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergyStar1 = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergyStar2 = new RDouble[numberOfSpecies];
    RDouble *transVibrationTime = new RDouble[numberOfSpecies];
    RDouble *transElectronTime = new RDouble[numberOfSpecies];
    RDouble **OmegaDerivative = NewPointer2 <RDouble> (numberOfSpecies, ntmodel - 1);

    //! Compute the partial derivatives in the main diagonal position.
    RDouble gasConstant = rjmk, mixGasCvtr, mixGasCvv, mixGasCve;
    RDouble Ttr = 0.0, Te = 0.0, electronPressure = 0.0, electronDensity = 0.0;
    RDouble vibrationOmega = 0.0, electronOmega = 0.0;
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermoProperty = molecularProperty->GetThermodynamicProperties();

    RDouble density = primitiveVariables[IR] * referenceDensity;
    RDouble pressure = primitiveVariables[IP] * referenceDensity * referenceVelocity * referenceVelocity;

    for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
    {
        //! Obtain the partial derivatives of chemical source to the vibration and electron energy.
        Omegas[s] = chemicalSourceTerms[s] / (referenceSpeciesEquation * cellVolume);
        for (int it = 1; it < ntmodel; ++ it)
        {
            OmegaDerivative[s][it - 1] = sourceDerivatives[s][numberOfSpecies + it - 1] / (referenceEnergyEquation * cellVolume);
        }
    }
    //! Obtain the mass fraction of species.
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        massFractions[s] = primitiveVariables[nm + s];
    }
    Ttr = temperatures[ITT] * referenceTemperature;
    Te  = temperatures[ITE] * referenceTemperature;

    //! Compute the specific heat at constant volume.
    mixGasCvtr = this->GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(massFractions);
    this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperatures[ITV], speciesCvv);
    this->ComputeDimensionalElectronSpecificHeatAtConstantVolume(temperatures[ITE], speciesCve);
    this->ComputeMixtureByMassFraction(massFractions, speciesCvv, mixGasCvv);
    this->ComputeMixtureByMassFraction(massFractions, speciesCve, mixGasCve);

    //! Compute vibration energy by translation-rotation temperature.
    this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperatures[ITT], speciesCvvStar1);
    //! Compute vibration energy by electron temperature.
    this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperatures[ITE], speciesCvvStar2);

    //! Compute vibration energy and electron energy.
    this->ComputeDimensionalVibrationEnergy(temperatures[ITV], vibrationEnergy);
    this->ComputeDimensionalElectronEnergy(temperatures[ITE], electronEnergy);
    this->ComputeDimensionalVibrationEnergy(temperatures[ITT], vibrationEnergyStar1);
    this->ComputeDimensionalVibrationEnergy(temperatures[ITE], vibrationEnergyStar2);

    //! Compute the relaxation time.of translation-vibration mode.
    this->GetTranslationAndVibrationRelaxationTime(thermoProperty, massFractions, speciesMass, Ttr, pressure, density, transVibrationTime);
    if (this->nElectronIndex >= 0)
    {
        electronDensity = density * massFractions[this->nElectronIndex];
        electronPressure = electronDensity * gasConstant * speciesMassReciprocal[this->nElectronIndex] * Te;
    }
    //! Compute the relaxation time.of electron-vibration mode.
    this->GetElectronAndVibrationRelaxationTime(thermoProperty, electronPressure, Te, transElectronTime);

    RDouble partv1 = 0.0, partv2 = 0.0, partv3 = 0.0, partv4 = 0.0;
    RDouble parte1 = 0.0, parte2 = 0.0, parte3 = 0.0;
    vibrationOmega = 0.0;
    electronOmega  = 0.0;
    for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
    {
        parte1 += Omegas[s] * speciesCve[s];
        parte2 += OmegaDerivative[s][1] * electronEnergy[s];
        //! Compute source term of electron energy equation.
        electronOmega += Omegas[s] * electronEnergy[s];

        if (thermoProperty[s].nType < 2)
        {
            continue;
        }
        partv1 += massFractions[s] * speciesCvvStar1[s] / transVibrationTime[s];
        partv2 += massFractions[s] * speciesCvv[s] / transVibrationTime[s];
        partv3 += Omegas[s] * speciesCvv[s];
        partv4 += OmegaDerivative[s][0] * vibrationEnergy[s];

        //! Compute source term of vibration energy equation.
        vibrationOmega += density * massFractions[s] * (vibrationEnergyStar1[s] - vibrationEnergy[s]) / transVibrationTime[s];
        vibrationOmega += Omegas[s] * vibrationEnergy[s];
    }
    if (nNitrogenIndex >= 0)
    {
        partv2 += massFractions[nNitrogenIndex] * speciesCvv[nNitrogenIndex] / transElectronTime[nNitrogenIndex];
        parte3 += massFractions[nNitrogenIndex] * speciesCvvStar2[nNitrogenIndex] / transElectronTime[nNitrogenIndex];
        //! Compute source term.
        RDouble temp = density * massFractions[nNitrogenIndex] * (vibrationEnergyStar2[nNitrogenIndex] - vibrationEnergy[nNitrogenIndex]) / transElectronTime[nNitrogenIndex];
        vibrationOmega += temp;
        electronOmega -= temp;
    }

    partv1 = partv1 / mixGasCvtr;
    partv2 = partv2 / mixGasCvv;
    partv3 = partv3 / (density * mixGasCvv);
    parte1 = parte1 / (density * mixGasCve);
    parte3 = parte3 / mixGasCve;

    //! Compute derivatives.
    sourceDerivatives[numberOfSpecies][numberOfSpecies] = partv3 + partv4 - partv1 - partv2;
    sourceDerivatives[numberOfSpecies + 1][numberOfSpecies + 1] = parte1 + parte2 - parte3;

    if (this->nElectronIndex >= 0)     //! Species including electron.
    {
        //! Compute the number density of electron.
        RDouble Ne = AVOGADRO_CONSTANT * electronDensity * speciesMassReciprocal[nElectronIndex];
        this->GetEffectiveCollisionCrossSectionalArea(thermoProperty, Te, Ne, sigma);

        RDouble Fe = 0.0, Ns = 0.0, fes = 0.0, speciesDensity = 0.0;
        RDouble heatVelocity = sqrt(8.0 * gasConstant * Te * speciesMassReciprocal[nElectronIndex] / PI);
        for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
        {
            speciesDensity = density * massFractions[s];
            Ns = AVOGADRO_CONSTANT * speciesDensity * speciesMassReciprocal[s];
            fes = Ns * sigma[s] * heatVelocity;
            Fe += fes * speciesMassReciprocal[s];
        }
        //! Compute source term of electron energy equation.
        electronOmega += 3.0 * electronDensity * gasConstant * (Ttr - Te) * Fe;

        RDouble mixGasEe = 0.0;
        this->ComputeMixtureByMassFraction(massFractions, electronEnergy , mixGasEe);
        sourceDerivatives[numberOfSpecies + 1][numberOfSpecies + 1] -= 3.0 * massFractions[nElectronIndex] * gasConstant * Fe * (1.0 / mixGasCvtr + 1.0 / mixGasCve + (Te - Ttr) / (2.0 * mixGasEe));
    }
    sourceDerivatives[numberOfSpecies][numberOfSpecies] *= reference1 * cellVolume;
    sourceDerivatives[numberOfSpecies + 1][numberOfSpecies + 1] *= reference1 * cellVolume;

    //! The source terms of vibration energy equation and electron energy equation.
    chemicalSourceTerms[numberOfSpecies] = reference2 * vibrationOmega * cellVolume;
    chemicalSourceTerms[numberOfSpecies + 1] = reference2 * electronOmega * cellVolume;

    delete []Omegas;
    delete []sigma;
    delete []speciesCvv;
    delete []speciesCve;
    delete []speciesCvvStar1;
    delete []speciesCvvStar2;
    delete []vibrationEnergy;
    delete []electronEnergy;
    delete []vibrationEnergyStar1;
    delete []vibrationEnergyStar2;
    delete []transVibrationTime;
    delete []transElectronTime;
    DelPointer2(OmegaDerivative);
}

void MixedGas::ComputeEnergySourceAndJacobianT2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{
    using namespace IDX;    //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble squareVelocity = referenceVelocity * referenceVelocity;
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble referenceEnergyEquation = referenceVelocity * reynoldsReferenceLengthDimensional;
    RDouble reference1 = reynoldsReferenceLengthDimensional / referenceVelocity;
    RDouble reference2 = referenceSpeciesEquation / squareVelocity;

    RDouble *Omegas = new RDouble[numberOfSpecies];
    RDouble *sigma = new RDouble[numberOfSpecies];
    RDouble *speciesCvv = new RDouble[numberOfSpecies];
    RDouble *speciesCve = new RDouble[numberOfSpecies];
    RDouble *speciesCvvStar1 = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergy = new RDouble[numberOfSpecies];
    RDouble *electronEnergy = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergyStar1 = new RDouble[numberOfSpecies];
    RDouble *transVibrationTime = new RDouble[numberOfSpecies];
    RDouble **OmegaDerivative = NewPointer2 <RDouble> (numberOfSpecies, ntmodel - 1);

    //! Compute the partial derivatives in the main diagonal position.
    RDouble mixGasCvtr, mixGasCvv, mixGasCve, mixGasCvve;
    RDouble Ttr = 0.0, Tve = 0.0, vibrationOmega = 0.0, temperature[3];
    RDouble gasConstant = rjmk, density = 0.0, pressure = 0.0;
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermoProperty = molecularProperty->GetThermodynamicProperties();

    for (int iCell = 0; iCell < nCellNumber; ++ iCell)
    {

        density = primitiveVariables[IR][iCell] * referenceDensity;
        pressure = primitiveVariables[IP][iCell] * referenceDensity * squareVelocity;

        for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
        {
            //! Obtain the partial derivatives of chemical source to the vibration and electron energy.
            Omegas[s] = chemicalSourceTerms[s][iCell] / (referenceSpeciesEquation * cellVolume[iCell]);
            for (int it = 1; it < ntmodel; ++ it)
            {
                OmegaDerivative[s][it - 1] = sourceDerivatives[s][numberOfSpecies + it - 1][iCell] / (referenceEnergyEquation * cellVolume[iCell]);
            }
        }
        //! Obtain the mass fraction of species.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            massFractions[s] = primitiveVariables[nm + s][iCell];
        }
        for (int it = 0; it < 3; ++ it)
        {
            temperature[it] = temperatures[it][iCell];
        }
        Ttr = temperature[ITT] * referenceTemperature;
        Tve = temperature[ITV] * referenceTemperature;

        //! Compute the specific heat at constant volume.
        mixGasCvtr = this->GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(massFractions);
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITV], speciesCvv);
        this->ComputeDimensionalElectronSpecificHeatAtConstantVolume(temperature[ITE], speciesCve);
        this->ComputeMixtureByMassFraction(massFractions, speciesCvv, mixGasCvv);
        this->ComputeMixtureByMassFraction(massFractions, speciesCve, mixGasCve);
        mixGasCvve = mixGasCvv + mixGasCve;

        //! Compute vibration energy by translation-rotation temperature.
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITT], speciesCvvStar1);
        //! Compute vibration energy and electron energy.
        this->ComputeDimensionalVibrationEnergy(temperature[ITV], vibrationEnergy);
        this->ComputeDimensionalElectronEnergy(temperature[ITE], electronEnergy);
        this->ComputeDimensionalVibrationEnergy(temperature[ITT], vibrationEnergyStar1);

        //! Compute the relaxation time.of translation-vibration mode.
        this->GetTranslationAndVibrationRelaxationTime(thermoProperty, massFractions, speciesMass, Ttr, pressure, density, transVibrationTime);

        RDouble part1 = 0.0, part2 = 0.0, part3 = 0.0, part4 = 0.0;
        vibrationOmega = 0.0;
        for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
        {
            if (thermoProperty[s].nType < 2)
            {
                continue;
            }

            part1 += massFractions[s] * speciesCvvStar1[s] / transVibrationTime[s];
            part2 += massFractions[s] * speciesCvv[s] / transVibrationTime[s];
            part3 += Omegas[s] * (speciesCvv[s] + speciesCve[s]);
            part4 += OmegaDerivative[s][0] * (vibrationEnergy[s] + electronEnergy[s]);

            //! Compute the source term of vibration energy equation.
            vibrationOmega += density * massFractions[s] * (vibrationEnergyStar1[s] - vibrationEnergy[s]) / transVibrationTime[s];
            vibrationOmega += Omegas[s] * (vibrationEnergy[s] + electronEnergy[s]);
        }
        part1 = part1 / mixGasCvtr;
        part2 = part2 / mixGasCvve;
        part3 = part3 / (density * mixGasCvve);
        //! Compute derivatives.
        RDouble vibrationJacobian = - part1 - part2 + part3 + part4;

        if (this->nElectronIndex >= 0) //species including electron.
        {
            //! Compute the number density of electron.
            RDouble electronDensity = density * massFractions[nElectronIndex];
            RDouble Ne = AVOGADRO_CONSTANT * electronDensity * speciesMassReciprocal[nElectronIndex];
            this->GetEffectiveCollisionCrossSectionalArea(thermoProperty, Tve, Ne, sigma);

            RDouble Fe = 0.0, Ns = 0.0, fes = 0.0, speciesDensity = 0.0;
            //! The heat motion velocity of electron.
            RDouble heatVelocity = sqrt(8.0 * gasConstant * Tve * speciesMassReciprocal[nElectronIndex] / PI);
            for (int s = 0; s <= nLastSpeciesIndex; ++ s)
            {
                speciesDensity = density * massFractions[s];
                Ns = AVOGADRO_CONSTANT * speciesDensity * speciesMassReciprocal[s];
                fes = Ns * sigma[s] * heatVelocity;
                Fe += fes * speciesMassReciprocal[s];
            }
            //! Compute the source term of vibration energy equation.
            vibrationOmega += 3.0 * electronDensity * Fe * gasConstant * (Ttr - Tve);

            if (nTEnergyModel >= 1)
            {
                vibrationJacobian += (1.5 * massFractions[nElectronIndex] * gasConstant * Fe / mixGasCvve) * (Ttr / Tve - 1.0);
                vibrationJacobian -= 3.0 * massFractions[nElectronIndex] * gasConstant * Fe * (1.0 / mixGasCvtr + 1.0 / mixGasCvve);
            }
            else
            {
                RDouble mixGasEv = 0.0, mixGasEe = 0.0, mixGasEve;
                this->ComputeMixtureByMassFraction(massFractions, vibrationEnergy, mixGasEv);
                this->ComputeMixtureByMassFraction(massFractions, electronEnergy , mixGasEe);
                mixGasEve = mixGasEv + mixGasEe;
                vibrationJacobian -= 3.0 * massFractions[nElectronIndex] * gasConstant * Fe * (1.0 / mixGasCvtr + 1.0 / mixGasCvve + (Tve - Ttr) / (2.0 * mixGasEve));
            }
        }
        sourceDerivatives[numberOfSpecies][numberOfSpecies][iCell] = reference1 * vibrationJacobian * cellVolume[iCell];
        //! The source term of vibration energy equation.
        chemicalSourceTerms[numberOfSpecies][iCell] = reference2 * vibrationOmega * cellVolume[iCell];

        // Inside of the loop.
    }

    delete []Omegas;
    delete []sigma;
    delete []speciesCvv;
    delete []speciesCve;
    delete []speciesCvvStar1;
    delete []vibrationEnergy;
    delete []electronEnergy;
    delete []vibrationEnergyStar1;
    delete []transVibrationTime;
    DelPointer2(OmegaDerivative);
}

void MixedGas::ComputeEnergySourceAndJacobianT2D3(int *Inzone, RDouble4D &primitiveVariables, RDouble4D &temperatures, RDouble4D &speciesCvvs, RDouble4D &speciesCves, RDouble4D &speciesEvs, RDouble4D &speciesEes, RDouble3D &totalCvtr, RDouble3D &totalCvv, RDouble3D &totalCve, RDouble3D &cellVolume, RDouble4D &chemicalSourceTerms, RDouble5D &sourceDerivatives)
{                                              
    using namespace IDX;    //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble squareVelocity = referenceVelocity * referenceVelocity;
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble referenceEnergyEquation = referenceVelocity * reynoldsReferenceLengthDimensional;
    RDouble reference1 = reynoldsReferenceLengthDimensional / referenceVelocity;
    RDouble reference2 = referenceSpeciesEquation / squareVelocity;

    RDouble *Omegas = Thermo_Energy_temparay.Cps;;
    RDouble *sigma = Thermo_Energy_temparay.Cvs;
    RDouble *speciesCvv = Thermo_Energy_temparay.Cvvs;
    RDouble *speciesCve = Thermo_Energy_temparay.Cves;
    RDouble *speciesCvvStar1 = Thermo_Energy_temparay.Cvtrs;
    RDouble *vibrationEnergy = Thermo_Energy_temparay.Etrs;
    RDouble *electronEnergy = Thermo_Energy_temparay.Evs;
    RDouble *vibrationEnergyStar1 = Thermo_Energy_temparay.Ees;
    RDouble *transVibrationTime = Thermo_Energy_temparay.Es;
    RDouble **OmegaDerivative = Thermo_Energy_temparay.OmegaDerivative;

    //! Compute the partial derivatives in the main diagonal position.
    RDouble mixGasCvtr, mixGasCvv, mixGasCve, mixGasCvve;
    RDouble Ttr = 0.0, Tve = 0.0, vibrationOmega = 0.0, temperature[3];
    RDouble gasConstant = rjmk, density = 0.0, pressure = 0.0;
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermoProperty = molecularProperty->GetThermodynamicProperties();

    int ist, ied, jst, jed, kst, ked;

    ist = Inzone[0];
    ied = Inzone[1];
    jst = Inzone[2];
    jed = Inzone[3];
    kst = Inzone[4];
    ked = Inzone[5];

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                density = primitiveVariables(i, j, k, IR) * referenceDensity;
                pressure = primitiveVariables(i, j, k, IP) * referenceDensity * squareVelocity;

                for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
                {
                    //! Obtain the partial derivatives of chemical source to the vibration and electron energy.
                    Omegas[s] = chemicalSourceTerms(i, j, k, s) / (referenceSpeciesEquation * cellVolume(i, j, k));
                    for (int it = 1; it < ntmodel; ++ it)
                    {
                        OmegaDerivative[s][it - 1] = sourceDerivatives(i, j, k, s, numberOfSpecies + it - 1) / ( referenceEnergyEquation * cellVolume(i, j, k) );
                    }
                }
                //! Obtain the mass fraction of species.
                for (int s = 0; s < numberOfSpecies; ++ s)
                {
                    massFractions[s] = primitiveVariables(i, j, k, nm + s);
                    speciesCvv[s] = speciesCvvs(i, j, k, s) / referenceSpecificHeat;
                    speciesCve[s] = speciesCves(i, j, k, s) / referenceSpecificHeat;

                    vibrationEnergy[s] = speciesEvs(i, j, k, s) * referenceEnergy;
                    electronEnergy[s]  = speciesEes(i, j, k, s) * referenceEnergy;
                }
                mixGasCvtr = totalCvtr(i, j, k) / referenceSpecificHeat;
                mixGasCvv  = totalCvv(i, j, k)  / referenceSpecificHeat;
                mixGasCve  = totalCve(i, j, k)  / referenceSpecificHeat;

                temperature[ITT] = temperatures(i, j, k, ITT);
                temperature[ITV] = temperatures(i, j, k, ITV);
                temperature[ITE] = temperature[ITV];

                Ttr = temperature[ITT] * referenceTemperature;
                Tve = temperature[ITV] * referenceTemperature; 

                mixGasCvve = mixGasCvv + mixGasCve;    

                //! Compute vibration energy by translation-rotation temperature.
                //this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITT], speciesCvvStar1);
                //this->ComputeDimensionalVibrationEnergy(temperature[ITT], vibrationEnergyStar1);

                ComputeCvvsAndEvsForSpecies(temperature[ITT] * referenceTemperature, speciesCvvStar1, vibrationEnergyStar1, 1);

                //! Compute the relaxation time.of translation-vibration mode.
                //this->GetTranslationAndVibrationRelaxationTime(thermoProperty, massFractions, speciesMass, Ttr, pressure, density, transVibrationTime);
                GetTranslationAndVibrationRelaxationTimeR(thermoProperty, massFractions, speciesMass, Ttr, pressure, density, transVibrationTime);

                RDouble part1 = 0.0, part2 = 0.0, part3 = 0.0, part4 = 0.0, part5 = 0.0, part6 = 0.0;
                vibrationOmega = 0.0;
                for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
                {
                    if (thermoProperty[s].nType < 2)
                    {
                        continue;
                    }

                    part1 += massFractions[s] * speciesCvvStar1[s] / transVibrationTime[s];
                    part2 += massFractions[s] * speciesCvv[s] / transVibrationTime[s];
                    part3 += Omegas[s] * (speciesCvv[s] + speciesCve[s]);
                    part4 += OmegaDerivative[s][0] * (vibrationEnergy[s] + electronEnergy[s]);

                    //! Compute the source term of vibration energy equation.
                    vibrationOmega += density * massFractions[s] * (vibrationEnergyStar1[s] - vibrationEnergy[s]) / transVibrationTime[s];
                    vibrationOmega += Omegas[s] * (vibrationEnergy[s] + electronEnergy[s]);
                }
                part1 = part1 / mixGasCvtr;
                part2 = part2 / mixGasCvve;
                part3 = part3 / (density * mixGasCvve);
                //! Compute derivatives.
               

                if (this->nElectronIndex >= 0) //species including electron.
                {
                    //! Compute the number density of electron.
                    RDouble electronDensity = density * massFractions[nElectronIndex];
                    RDouble Ne = AVOGADRO_CONSTANT * electronDensity * speciesMassReciprocal[nElectronIndex];
                    this->GetEffectiveCollisionCrossSectionalArea(thermoProperty, Tve, Ne, sigma);

                    RDouble Fe = 0.0, Ns = 0.0, fes = 0.0, speciesDensity = 0.0;
                    //! The heat motion velocity of electron.
                    RDouble heatVelocity = sqrt(8.0 * gasConstant * Tve * speciesMassReciprocal[nElectronIndex] / PI);
                    for (int s = 0; s <= nLastSpeciesIndex; ++ s)
                    {
                        speciesDensity = density * massFractions[s];
                        Ns = AVOGADRO_CONSTANT * speciesDensity * speciesMassReciprocal[s];
                        fes = Ns * sigma[s] * heatVelocity;
                        Fe += fes * speciesMassReciprocal[s];
                    }
                    //! Compute the source term of vibration energy equation.
                    vibrationOmega += 3.0 * electronDensity * Fe * gasConstant * (Ttr - Tve);

                    if (nTEnergyModel >= 1 || nChemcalSourceModified > 0)
                    {
                        part5 = (1.5 * massFractions[nElectronIndex] * gasConstant * Fe / mixGasCvve) * (Ttr / Tve - 1.0);
                        part6 = 3.0 * massFractions[nElectronIndex] * gasConstant * Fe * (1.0 / mixGasCvtr + 1.0 / mixGasCvve);
                    }
                    else
                    {
                        RDouble mixGasEv = 0.0, mixGasEe = 0.0, mixGasEve;
                        this->ComputeMixtureByMassFraction(massFractions, vibrationEnergy, mixGasEv);
                        this->ComputeMixtureByMassFraction(massFractions, electronEnergy , mixGasEe);
                        mixGasEve = mixGasEv + mixGasEe;
                        //vibrationJacobian -= 3.0 * massFractions[nElectronIndex] * gasConstant * Fe * (1.0 / mixGasCvtr + 1.0 / mixGasCvve + (Tve - Ttr) / (2.0 * mixGasEve));
                        //part5 = 0.0;
                        part6 = 3.0 * massFractions[nElectronIndex] * gasConstant * Fe * (1.0 / mixGasCvtr + 1.0 / mixGasCvve + (Tve - Ttr) / (2.0 * mixGasEve));
                    }
                }

                 RDouble vibrationJacobian;
                if(nChemcalSourceModified == 0)
                {
                    vibrationJacobian = - part1 - part2 - part6 + part3 + part4 + part5;
                }
                else if(nChemcalSourceModified == 1)
                {
                    vibrationJacobian = - part1 - part2 - part6 + min(0.0, part3) + min(0.0, part4) + min(0.0, part5);
                }
                else
                {
                    vibrationJacobian = - abs(part1) - abs(part2) - abs(part6) -abs(part3) - abs(part4) - abs(part5);
                }
                sourceDerivatives(i, j, k, numberOfSpecies, numberOfSpecies) = reference1 * vibrationJacobian * cellVolume(i, j, k);
                //! The source term of vibration energy equation.
                chemicalSourceTerms(i, j, k, numberOfSpecies) = reference2 * vibrationOmega * cellVolume(i, j, k);
            }
        }
        // Inside of the loop.
    }
}


void MixedGas::ComputeEnergySourceAndJacobianT3(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{
    using namespace IDX;     //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble squareVelocity = referenceVelocity * referenceVelocity;
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble referenceEnergyEquation = referenceVelocity * reynoldsReferenceLengthDimensional;
    RDouble reference1 = reynoldsReferenceLengthDimensional / referenceVelocity;
    RDouble reference2 = referenceSpeciesEquation / squareVelocity;

    RDouble *Omegas = new RDouble[numberOfSpecies];
    RDouble *sigma = new RDouble[numberOfSpecies];
    RDouble *speciesCvv = new RDouble[numberOfSpecies];
    RDouble *speciesCve = new RDouble[numberOfSpecies];
    RDouble *speciesCvvStar1 = new RDouble[numberOfSpecies];
    RDouble *speciesCvvStar2 = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergy = new RDouble[numberOfSpecies];
    RDouble *electronEnergy = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergyStar1 = new RDouble[numberOfSpecies];
    RDouble *vibrationEnergyStar2 = new RDouble[numberOfSpecies];
    RDouble *transVibrationTime = new RDouble[numberOfSpecies];
    RDouble *transElectronTime = new RDouble[numberOfSpecies];
    RDouble **OmegaDerivative = NewPointer2 <RDouble> (numberOfSpecies, ntmodel - 1);

    //! Compute the partial derivatives in the main diagonal position.
    RDouble gasConstant = rjmk, mixGasCvtr, mixGasCvv, mixGasCve;
    RDouble Ttr = 0.0, Te = 0.0, electronPressure = 0.0, electronDensity = 0.0;
    RDouble vibrationOmega = 0.0, electronOmega = 0.0, temperature[3];
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermoProperty = molecularProperty->GetThermodynamicProperties();

    for (int iCell = 0; iCell < nCellNumber; ++ iCell)
    {

        RDouble density = primitiveVariables[IR][iCell] * referenceDensity;
        RDouble pressure = primitiveVariables[IP][iCell] * referenceDensity * squareVelocity;

        for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
        {
            //! Obtain the partial derivatives of chemical source to the vibration and electron energy.
            Omegas[s] = chemicalSourceTerms[s][iCell] / (referenceSpeciesEquation * cellVolume[iCell]);
            for (int it = 1; it < ntmodel; ++ it)
            {
                OmegaDerivative[s][it - 1] = sourceDerivatives[s][numberOfSpecies + it - 1][iCell] / (referenceEnergyEquation * cellVolume[iCell]);
            }
        }
        //! Obtain the mass fraction of species.
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            massFractions[s] = primitiveVariables[nm + s][iCell];
        }
        for (int it = 0; it < 3; ++ it)
        {
            temperature[it] = temperatures[it][iCell];
        }
        Ttr = temperature[ITT] * referenceTemperature;
        Te  = temperature[ITE] * referenceTemperature;

        //! Compute the specific heat at constant volume.
        mixGasCvtr = this->GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(massFractions);
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITV], speciesCvv);
        this->ComputeDimensionalElectronSpecificHeatAtConstantVolume(temperature[ITE], speciesCve);
        this->ComputeMixtureByMassFraction(massFractions, speciesCvv, mixGasCvv);
        this->ComputeMixtureByMassFraction(massFractions, speciesCve, mixGasCve);

        //! Compute vibration energy by translation-rotation temperature.
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITT], speciesCvvStar1);
        //! Compute vibration energy by electron temperature.
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITE], speciesCvvStar2);

        //! Compute vibration energy and electron energy.
        this->ComputeDimensionalVibrationEnergy(temperature[ITV], vibrationEnergy);
        this->ComputeDimensionalElectronEnergy(temperature[ITE], electronEnergy);
        this->ComputeDimensionalVibrationEnergy(temperature[ITT], vibrationEnergyStar1);
        this->ComputeDimensionalVibrationEnergy(temperature[ITE], vibrationEnergyStar2);

        //! Compute the relaxation time.of translation-vibration mode.
        this->GetTranslationAndVibrationRelaxationTime(thermoProperty, massFractions, speciesMass, Ttr, pressure, density, transVibrationTime);
        if (this->nElectronIndex >= 0)
        {
            electronDensity = density * massFractions[this->nElectronIndex];
            electronPressure = electronDensity * gasConstant * speciesMassReciprocal[this->nElectronIndex] * Te;
        }
        //! Compute the relaxation time.of electron-vibration mode.
        this->GetElectronAndVibrationRelaxationTime(thermoProperty, electronPressure, Te, transElectronTime);

        RDouble partv1 = 0.0, partv2 = 0.0, partv3 = 0.0, partv4 = 0.0;
        RDouble parte1 = 0.0, parte2 = 0.0, parte3 = 0.0;
        vibrationOmega = 0.0;
        electronOmega  = 0.0;
        for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
        {
            parte1 += Omegas[s] * speciesCve[s];
            parte2 += OmegaDerivative[s][1] * electronEnergy[s];
            //! Compute source term of electron energy equation.
            electronOmega += Omegas[s] * electronEnergy[s];

            if (thermoProperty[s].nType < 2)
            {
                continue;
            }
            partv1 += massFractions[s] * speciesCvvStar1[s] / transVibrationTime[s];
            partv2 += massFractions[s] * speciesCvv[s] / transVibrationTime[s];
            partv3 += Omegas[s] * speciesCvv[s];
            partv4 += OmegaDerivative[s][0] * vibrationEnergy[s];

            //! Compute source term of vibration energy equation.
            vibrationOmega += density * massFractions[s] * (vibrationEnergyStar1[s] - vibrationEnergy[s]) / transVibrationTime[s];
            vibrationOmega += Omegas[s] * vibrationEnergy[s];
        }
        if (nNitrogenIndex >= 0)
        {
            partv2 += massFractions[nNitrogenIndex] * speciesCvv[nNitrogenIndex] / transElectronTime[nNitrogenIndex];
            parte3 += massFractions[nNitrogenIndex] * speciesCvvStar2[nNitrogenIndex] / transElectronTime[nNitrogenIndex];
            //! Compute source term.
            RDouble temp = density * massFractions[nNitrogenIndex] * (vibrationEnergyStar2[nNitrogenIndex] - vibrationEnergy[nNitrogenIndex]) / transElectronTime[nNitrogenIndex];
            vibrationOmega += temp;
            electronOmega -= temp;
        }

        partv1 = partv1 / mixGasCvtr;
        partv2 = partv2 / mixGasCvv;
        partv3 = partv3 / (density * mixGasCvv);
        parte1 = parte1 / (density * mixGasCve);
        parte3 = parte3 / mixGasCve;

        //! Compute derivatives.
        RDouble vibrationJacobian = partv3 + partv4 - partv1 - partv2;
        RDouble electronJacobian  = parte1 + parte2 - parte3;

        if (this->nElectronIndex >= 0)     //! Species including electron.
        {
            //! Compute the number density of electron.
            RDouble Ne = AVOGADRO_CONSTANT * electronDensity * speciesMassReciprocal[nElectronIndex];
            this->GetEffectiveCollisionCrossSectionalArea(thermoProperty, Te, Ne, sigma);

            RDouble Fe = 0.0, Ns = 0.0, fes = 0.0, speciesDensity = 0.0;
            RDouble heatVelocity = sqrt(8.0 * gasConstant * Te * speciesMassReciprocal[nElectronIndex] / PI);
            for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
            {
                speciesDensity = density * massFractions[s];
                Ns = AVOGADRO_CONSTANT * speciesDensity * speciesMassReciprocal[s];
                fes = Ns * sigma[s] * heatVelocity;
                Fe += fes * speciesMassReciprocal[s];
            }
            //! Compute source term of electron energy equation.
            electronOmega += 3.0 * electronDensity * Fe * gasConstant * (Ttr - Te) * Fe;

            RDouble mixGasEe = 0.0;
            this->ComputeMixtureByMassFraction(massFractions, electronEnergy , mixGasEe);
            electronJacobian -= 3.0 * massFractions[nElectronIndex] * gasConstant * Fe * (1.0 / mixGasCvtr + 1.0 / mixGasCve + (Te - Ttr) / (2.0 * mixGasEe));
        }
        sourceDerivatives[numberOfSpecies][numberOfSpecies][iCell] = reference1 * vibrationJacobian * cellVolume[iCell];
        sourceDerivatives[numberOfSpecies + 1][numberOfSpecies + 1][iCell] = reference1 * electronJacobian * cellVolume[iCell];

        //! The source terms of vibration energy equation and electron energy equation.
        chemicalSourceTerms[numberOfSpecies][iCell] = reference2 * vibrationOmega * cellVolume[iCell];
        chemicalSourceTerms[numberOfSpecies + 1][iCell] = reference2 * electronOmega * cellVolume[iCell];

        //Inside of the loop.
    }

    delete []Omegas;
    delete []sigma;
    delete []speciesCvv;
    delete []speciesCve;
    delete []speciesCvvStar1;
    delete []speciesCvvStar2;
    delete []vibrationEnergy;
    delete []electronEnergy;
    delete []vibrationEnergyStar1;
    delete []vibrationEnergyStar2;
    delete []transVibrationTime;
    delete []transElectronTime;
    DelPointer2(OmegaDerivative);
}

void MixedGas::ComputeEnergySourceAndJacobianT3D3(int *Inzone, RDouble4D &primitiveVariables, RDouble4D &temperatures, RDouble4D &speciesCvvs, RDouble4D &speciesCves, RDouble4D &speciesEvs, RDouble4D &speciesEes, RDouble3D &totalCvtr, RDouble3D &totalCvv, RDouble3D &totalCve, RDouble3D &cellVolume, RDouble4D &chemicalSourceTerms, RDouble5D &sourceDerivatives)
{
    using namespace IDX;     //! The following variables is dimensional values.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocalDimensional();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceDensity = this->referenceParameterDimensional.GetDensity();
    RDouble squareVelocity = referenceVelocity * referenceVelocity;
    RDouble referenceSpeciesEquation = reynoldsReferenceLengthDimensional / (referenceDensity * referenceVelocity);
    RDouble referenceEnergyEquation = referenceVelocity * reynoldsReferenceLengthDimensional;
    RDouble reference1 = reynoldsReferenceLengthDimensional / referenceVelocity;
    RDouble reference2 = referenceSpeciesEquation / squareVelocity;

    //RDouble *Omegas = new RDouble[numberOfSpecies];
    //RDouble *sigma = new RDouble[numberOfSpecies];
    //RDouble *speciesCvv = new RDouble[numberOfSpecies];
    //RDouble *speciesCve = new RDouble[numberOfSpecies];
    //RDouble *speciesCvvStar1 = new RDouble[numberOfSpecies];
    //RDouble *speciesCvvStar2 = new RDouble[numberOfSpecies];
    //RDouble *vibrationEnergy = new RDouble[numberOfSpecies];
    //RDouble *electronEnergy = new RDouble[numberOfSpecies];
    //RDouble *vibrationEnergyStar1 = new RDouble[numberOfSpecies];
    //RDouble *vibrationEnergyStar2 = new RDouble[numberOfSpecies];
    //RDouble *transVibrationTime = new RDouble[numberOfSpecies];
    //RDouble *transElectronTime = new RDouble[numberOfSpecies];
    //RDouble **OmegaDerivative = NewPointer2 <RDouble> (numberOfSpecies, ntmodel - 1);

    RDouble *Omegas = Thermo_Energy_temparay.Cps;;
    RDouble *sigma = Thermo_Energy_temparay.Cvs;
    RDouble *speciesCvv = Thermo_Energy_temparay.Cvvs;
    RDouble *speciesCve = Thermo_Energy_temparay.Cves;
    RDouble *speciesCvvStar1 = Thermo_Energy_temparay.Cvtrs;
    RDouble *speciesCvvStar2 = Thermo_Energy_temparay.Rms;
    RDouble *vibrationEnergy = Thermo_Energy_temparay.Etrs;
    RDouble *electronEnergy = Thermo_Energy_temparay.Evs;
    RDouble *vibrationEnergyStar1 = Thermo_Energy_temparay.Ees;
    RDouble *vibrationEnergyStar2 = Thermo_Energy_temparay.Hs;
    RDouble *transVibrationTime = Thermo_Energy_temparay.Es;
    RDouble *transElectronTime = workSpecies;
    RDouble **OmegaDerivative = Thermo_Energy_temparay.OmegaDerivative;

    //! Compute the partial derivatives in the main diagonal position.
    RDouble gasConstant = rjmk, mixGasCvtr, mixGasCvv, mixGasCve;
    RDouble Ttr = 0.0, Te = 0.0, electronPressure = 0.0, electronDensity = 0.0;
    RDouble vibrationOmega = 0.0, electronOmega = 0.0, temperature[3];
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    Thermo_Param *thermoProperty = molecularProperty->GetThermodynamicProperties();

    int ist, ied, jst, jed, kst, ked;
    ist = Inzone[0];
    ied = Inzone[1];
    jst = Inzone[2];
    jed = Inzone[3];
    kst = Inzone[4];
    ked = Inzone[5];

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble density = primitiveVariables(i, j, k, IR) * referenceDensity;
                RDouble pressure = primitiveVariables(i, j, k, IP) * referenceDensity * squareVelocity;

                for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
                {
                    //! Obtain the partial derivatives of chemical source to the vibration and electron energy.
                    Omegas[s] = chemicalSourceTerms(i, j, k, s) / (referenceSpeciesEquation * cellVolume(i, j, k));
                    for (int it = 1; it < ntmodel; ++ it)
                    {
                        OmegaDerivative[s][it - 1] = sourceDerivatives(i, j, k, s, numberOfSpecies + it - 1) / (referenceEnergyEquation * cellVolume(i, j, k));
                    }
                }
                //! Obtain the mass fraction of species.
                for (int s = 0; s < numberOfSpecies; ++ s)
                {
                    massFractions[s] = primitiveVariables(i, j, k, nm + s);
                    speciesCvv[s] = speciesCvvs(i, j, k, s) / referenceSpecificHeat;
                    speciesCve[s] = speciesCves(i, j, k, s) / referenceSpecificHeat;

                    vibrationEnergy[s] = speciesEvs(i, j, k, s) * referenceEnergy;
                    electronEnergy[s]  = speciesEes(i, j, k, s) * referenceEnergy;
                }
                mixGasCvtr = totalCvtr(i, j, k) / referenceSpecificHeat;
                mixGasCvv  = totalCvv(i, j, k)  / referenceSpecificHeat;
                mixGasCve  = totalCve(i, j, k)  / referenceSpecificHeat;

                for (int it = 0; it < 3; ++ it)
                {
                    temperature[it] = temperatures(i, j, k, it);
                }
                Ttr = temperature[ITT] * referenceTemperature;
                Te  = temperature[ITE] * referenceTemperature;

                //! Compute the specific heat at constant volume.
                /*mixGasCvtr = this->GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(massFractions);
                this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITV], speciesCvv);
                this->ComputeDimensionalElectronSpecificHeatAtConstantVolume(temperature[ITE], speciesCve);
                this->ComputeMixtureByMassFraction(massFractions, speciesCvv, mixGasCvv);
                this->ComputeMixtureByMassFraction(massFractions, speciesCve, mixGasCve);*/

                //! Compute vibration energy by translation-rotation temperature.
                this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITT], speciesCvvStar1);
                //! Compute vibration energy by electron temperature.
                this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(temperature[ITE], speciesCvvStar2);

                //! Compute vibration energy and electron energy.
                /*this->ComputeDimensionalVibrationEnergy(temperature[ITV], vibrationEnergy);
                this->ComputeDimensionalElectronEnergy(temperature[ITE], electronEnergy);*/
                this->ComputeDimensionalVibrationEnergy(temperature[ITT], vibrationEnergyStar1);
                this->ComputeDimensionalVibrationEnergy(temperature[ITE], vibrationEnergyStar2);

                //! Compute the relaxation time.of translation-vibration mode.
                this->GetTranslationAndVibrationRelaxationTime(thermoProperty, massFractions, speciesMass, Ttr, pressure, density, transVibrationTime);
                if (this->nElectronIndex >= 0)
                {
                    electronDensity = density * massFractions[this->nElectronIndex];
                    electronPressure = electronDensity * gasConstant * speciesMassReciprocal[this->nElectronIndex] * Te;
                }
                //! Compute the relaxation time.of electron-vibration mode.
                this->GetElectronAndVibrationRelaxationTime(thermoProperty, electronPressure, Te, transElectronTime);

                RDouble partv1 = 0.0, partv2 = 0.0, partv3 = 0.0, partv4 = 0.0;
                RDouble parte1 = 0.0, parte2 = 0.0, parte3 = 0.0;
                vibrationOmega = 0.0;
                electronOmega  = 0.0;
                for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
                {
                    parte1 += Omegas[s] * speciesCve[s];
                    parte2 += OmegaDerivative[s][1] * electronEnergy[s];
                    //! Compute source term of electron energy equation.
                    electronOmega += Omegas[s] * electronEnergy[s];

                    if (thermoProperty[s].nType < 2)
                    {
                        continue;
                    }
                    partv1 += massFractions[s] * speciesCvvStar1[s] / transVibrationTime[s];
                    partv2 += massFractions[s] * speciesCvv[s] / transVibrationTime[s];
                    partv3 += Omegas[s] * speciesCvv[s];
                    partv4 += OmegaDerivative[s][0] * vibrationEnergy[s];

                    //! Compute source term of vibration energy equation.
                    vibrationOmega += density * massFractions[s] * (vibrationEnergyStar1[s] - vibrationEnergy[s]) / transVibrationTime[s];
                    vibrationOmega += Omegas[s] * vibrationEnergy[s];
                }
                if (nNitrogenIndex >= 0)
                {
                    partv2 += massFractions[nNitrogenIndex] * speciesCvv[nNitrogenIndex] / transElectronTime[nNitrogenIndex];
                    parte3 += massFractions[nNitrogenIndex] * speciesCvvStar2[nNitrogenIndex] / transElectronTime[nNitrogenIndex];
                    //! Compute source term.
                    RDouble temp = density * massFractions[nNitrogenIndex] * (vibrationEnergyStar2[nNitrogenIndex] - vibrationEnergy[nNitrogenIndex]) / transElectronTime[nNitrogenIndex];
                    vibrationOmega += temp;
                    electronOmega -= temp;
                }

                partv1 = partv1 / mixGasCvtr;
                partv2 = partv2 / mixGasCvv;
                partv3 = partv3 / (density * mixGasCvv);
                parte1 = parte1 / (density * mixGasCve);
                parte3 = parte3 / mixGasCve;

                //! Compute derivatives.
                RDouble vibrationJacobian = partv3 + partv4 - partv1 - partv2;
                RDouble electronJacobian  = parte1 + parte2 - parte3;

                if (this->nElectronIndex >= 0)     //! Species including electron.
                {
                    //! Compute the number density of electron.
                    RDouble Ne = AVOGADRO_CONSTANT * electronDensity * speciesMassReciprocal[nElectronIndex];
                    this->GetEffectiveCollisionCrossSectionalArea(thermoProperty, Te, Ne, sigma);

                    RDouble Fe = 0.0, Ns = 0.0, fes = 0.0, speciesDensity = 0.0;
                    RDouble heatVelocity = sqrt(8.0 * gasConstant * Te * speciesMassReciprocal[nElectronIndex] / PI);
                    for (int s = 0; s <= this->nLastSpeciesIndex; ++ s)
                    {
                        speciesDensity = density * massFractions[s];
                        Ns = AVOGADRO_CONSTANT * speciesDensity * speciesMassReciprocal[s];
                        fes = Ns * sigma[s] * heatVelocity;
                        Fe += fes * speciesMassReciprocal[s];
                    }
                    //! Compute source term of electron energy equation.
                    electronOmega += 3.0 * electronDensity * Fe * gasConstant * (Ttr - Te) * Fe;

                    RDouble mixGasEe = 0.0;
                    this->ComputeMixtureByMassFraction(massFractions, electronEnergy , mixGasEe);
                    electronJacobian -= 3.0 * massFractions[nElectronIndex] * gasConstant * Fe * (1.0 / mixGasCvtr + 1.0 / mixGasCve + (Te - Ttr) / (2.0 * mixGasEe));
                }
                sourceDerivatives(i, j, k, numberOfSpecies, numberOfSpecies) = reference1 * vibrationJacobian * cellVolume(i, j, k);
                sourceDerivatives(i, j, k, numberOfSpecies + 1, numberOfSpecies + 1)  = reference1 * electronJacobian * cellVolume(i, j, k);

                //! The source terms of vibration energy equation and electron energy equation.
                chemicalSourceTerms(i, j, k, numberOfSpecies) = reference2 * vibrationOmega * cellVolume(i, j, k);
                chemicalSourceTerms(i, j, k, numberOfSpecies + 1) = reference2 * electronOmega * cellVolume(i, j, k);

                //Inside of the loop.
            }
        }
    }
}

void MixedGas::GetLinearlyAveragedPolynomialCoefficients(int nSpecies, RDouble temperature, RDouble *polynomialCoef)
{
    int nIndex1, nIndex2;      //! Record the index of temperature range.
    RDouble coef1, coef2;      //! The coefficients for linear polynomial on the temperature boundary.
    RDouble T = temperature;   //! temperature is the dimensional value.

    //! Select the temperature range and compute the coefficients.
    if (T > 15500.0)
    {
        nIndex1 = 3;
        nIndex2 = 4;
        if (T > 25500.0)     //! Range: (25500, 30000).
        {
            coef1 = 0.0;
            coef2 = 1.0;
        }
        else if (T <= 24500.0)     //! Range: (15500, 24500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else     //! Range: (24500, 25500).
        {
            coef2 = 0.001 * (T - 24500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else if (T > 6500.0)
    {
        nIndex1 = 2;
        nIndex2 = 3;
        if (T <= 14500.0)     //! Range: (6500, 14500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (14500, 15500).
        {
            coef2 = 0.001 * (T - 14500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else if (T > 1200.0)
    {
        nIndex1 = 1;
        nIndex2 = 2;
        if (T <= 5500.0)     //! Range: (1200, 5500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (5500, 6500).
        {
            coef2 = 0.001 * (T - 5500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else
    {
        nIndex1 = 0;
        nIndex2 = 1;
        if (T <= 800.0)      //! Range: (0, 800).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (800, 1200).
        {
            coef2 = (T - 800.0)/400.0;
            coef1 = 1.0 - coef2;
        }
    }

    RDouble *a = this->thermodynamicManager->GetPolynomialCoefficient(nSpecies, nIndex1);
    RDouble *b = this->thermodynamicManager->GetPolynomialCoefficient(nSpecies, nIndex2);
    int n = thermodynamicManager->GetNumberofPolynomialCoefficient();
    for (int i = 0; i < n; ++ i)
    {
        polynomialCoef[i] = coef1 * a[i] + coef2 * b[i];
    }
}

void MixedGas::GetLinearlyAveragedPolynomialCoefficientsT2(RDouble temperature, RDouble **polynomialCoef)
{
    int nIndex1, nIndex2;      //! Record the index of temperature range.
    RDouble coef1, coef2;      //! The coefficients for linear polynomial on the temperature boundary.
    RDouble T = temperature;   //! temperature is the dimensional value.

    //! Select the temperature range and compute the coefficients.
    if (T > 15500.0)
    {
        nIndex1 = 4;
        nIndex2 = 5;
        if (T > 25500.0)     //! Range: (25500, 30000).
        {
            coef1 = 0.0;
            coef2 = 1.0;
        }
        else if (T <= 24500.0)     //! Range: (15500, 24500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else     //! Range: (24500, 25500).
        {
            coef2 = 0.001 * (T - 24500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else if (T > 6500.0)
    {
        nIndex1 = 3;
        nIndex2 = 4;
        if (T <= 14500.0)     //! Range: (6500, 14500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (14500, 15500).
        {
            coef2 = 0.001 * (T - 14500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else if (T > 1200.0)
    {
        nIndex1 = 2;
        nIndex2 = 3;
        if (T <= 5500.0)     //! Range: (1200, 5500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (5500, 6500).
        {
            coef2 = 0.001 * (T - 5500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else
    {
        nIndex1 = 1;
        nIndex2 = 2;
        if (T <= 800.0)      //! Range: (0, 800).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (800, 1200).
        {
            coef2 = (T - 800.0)/400.0;
            coef1 = 1.0 - coef2;
        }
    }

    RDouble benchmarkTemp = 200, T1 = 0.0, T2 = 0.0;
    RDouble *polynomial1 = NULL, *polynomial2 = NULL;

    int n = thermodynamicManager->GetNumberofPolynomialCoefficient();
    if (ABS(coef1 - 1) <= EPSILON)
    {
        for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
        {
            polynomial1 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, nIndex1);
            for (int i = 0; i < n; ++ i)
            {
                polynomialCoef[iSpecies][i] = polynomial1[i];
            }
        }
    }
    else if (ABS(coef2 - 1) <= EPSILON)
    {
        for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
        {
            polynomial2 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, nIndex2);
            for (int i = 0; i < n; ++ i)
            {
                polynomialCoef[iSpecies][i] = polynomial2[i];
            }
        }
    }
    else
    {
        for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
        {
            polynomial1 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, nIndex1);
            polynomial2 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, nIndex2);
            for (int i = 0; i < n; ++ i)
            {
                polynomialCoef[iSpecies][i] = coef1 * polynomial1[i] + coef2 * polynomial2[i];
            }
        }
    }

    //Modify the coefficient when the temperature is low.
    bool isSmooth = false;
    for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
    {
        benchmarkTemp = this->thermodynamicManager->GetBenchmarkTemperature(iSpecies);
        if (isSmooth)
        {
            T1 = benchmarkTemp - 100;
            T2 = benchmarkTemp + 100;
            if (T < T1)
            {
                polynomial1 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, 0);
                for (int i = 0; i < n; ++ i)
                {
                    polynomialCoef[iSpecies][i] = polynomial1[i];
                }
            }
            else if (T <= T2)
            {
                polynomial1 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, 0);
                polynomial2 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, 1);
                coef2 = (T - T1) / 200.0;
                coef1 = 1.0 - coef2;
                for (int i = 0; i < n; ++ i)
                {
                    polynomialCoef[iSpecies][i] = coef1 * polynomial1[i] + coef2 * polynomial2[i];
                }
            }
        }
        else //No smooth.
        {
            if (T < benchmarkTemp)
            {
                polynomial1 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, 0);
                for (int i = 0; i < n; ++ i)
                {
                    polynomialCoef[iSpecies][i] = polynomial1[i];
                }
            }
        }
    }
}


void MixedGas::GetLinearlyAveragedPolynomialCoefficientsIndex(RDouble temperature,RDouble &coef, int &nIndex1, int &nIndex2)
{
    //int nIndex1, nIndex2;      //! Record the index of temperature range.
    RDouble coef1 = 1.0, coef2 = 0.0;      //! The coefficients for linear polynomial on the temperature boundary.
    RDouble T = temperature;   //! temperature is the dimensional value.

    //! Select the temperature range and compute the coefficients.
    if (T > 15500.0)
    {
        if (T > 25500.0)     //! Range: (25500, 30000).
        {
            nIndex1 = 4;
            nIndex2 = 0;
            coef = -1.0;
        }
        else if (T <= 24500.0)     //! Range: (15500, 24500).
        {
            nIndex1 = 3;
            nIndex2 = 0;
            coef = -1.0;
        }
        else     //! Range: (24500, 25500).
        {
            coef = 0.001 * (T - 24500.0);
            nIndex1 = 3;
            nIndex2 = 4;
        }
    }
    else if (T > 6500.0)
    {
        nIndex1 = 2;
        nIndex2 = 3;
        if (T <= 14500.0)     //! Range: (6500, 14500).
        {
            coef = -1.0;
        }
        else      //! Range: (14500, 15500).
        {
            coef = 0.001 * (T - 14500.0);
        }
    }
    else if (T > 1200.0)
    {
        nIndex1 = 1;
        nIndex2 = 2;
        if (T <= 5500.0)     //! Range: (1200, 5500).
        {
            coef = -1.0;
        }
        else      //! Range: (5500, 6500).
        {
            coef = 0.001 * (T - 5500.0);
        }
    }
    else
    {
        nIndex1 = 0;
        nIndex2 = 1;
        if (T <= 800.0)      //! Range: (0, 800).
        {
            coef = -1.0;
        }
        else      //! Range: (800, 1200).
        {
            coef = (T - 800.0)/400.0;
        }
    }
}


void MixedGas::GetLinearlyAveragedPolynomialCoefficientsByIndex(int nSpecies, RDouble *polynomialCoef, RDouble coef, int nIndex1, int nIndex2)
{
    RDouble *a = this->thermodynamicManager->GetPolynomialCoefficient(nSpecies, nIndex1);
    int n = thermodynamicManager->GetNumberofPolynomialCoefficient();
    if(coef < 0)
    {
        for (int i = 0; i < n; ++ i)
        {
            polynomialCoef[i] = a[i];
        }
    }
    else
    {
        RDouble *b = this->thermodynamicManager->GetPolynomialCoefficient(nSpecies, nIndex2);
        for (int i = 0; i < n; ++ i)
        {
            polynomialCoef[i] = (1.0 - coef) * a[i] + coef * b[i];
        }
    }
}
void MixedGas::GetAirInfo(const RDouble &h1, RDouble &t, RDouble &p, RDouble &den, RDouble &a)
{
    //! USSA 76 Model for atmosphere.
    RDouble h      = h1 * 1000.0;
    RDouble r      = 287.053;
    RDouble g0     = 9.80665;
    RDouble rp     = 6.37111e6;
    RDouble g      = (rp/(rp+h))*(rp/(rp+h))*g0;
    RDouble t0     = 288.15;
    RDouble p0     = 10.1325e2;
    RDouble rho0   = 1.225;
    RDouble t11    = 216.65;
    RDouble p11    = 2.2632e2;
    RDouble rho11  = 3.6392e-1;
    RDouble t20    = t11;
    RDouble p20    = 5.4747e1;
    RDouble rho20  = 8.8035e-2;
    RDouble t32    = 228.65;
    RDouble p32    = 8.6789;
    RDouble rho32  = 1.3225e-2;
    RDouble t47    = 270.65;
    RDouble p47    = 1.1090;
    RDouble rho47  = 1.4275e-3;
    RDouble t52    = t47;
    RDouble p52    = 5.8997e-1;
    RDouble rho52  = 7.5943e-4;
    RDouble t61    = 252.65;
    RDouble p61    = 1.8209e-1;
    RDouble rho61  = 2.5109e-4;
    RDouble t79    = 180.65;
    RDouble p79    = 1.0376e-2;
    RDouble rho79  = 2.0010e-5;
    RDouble t90    = t79;
    RDouble p90    = 1.6437e-3;
    RDouble rho90  = 3.4165e-6;
    RDouble t100   = 210.02;
    RDouble p100   = 3.0070e-4;
    RDouble rho100 = 5.6044e-7;
    RDouble t110   = 257.00;
    RDouble p110   = 7.3527e-5;
    RDouble rho110 = 9.7081e-8;
    RDouble t120   = 349.49;
    RDouble p120   = 2.5209e-5;
    RDouble rho120 = 2.2222e-8;
    RDouble t150   = 892.79;
    RDouble p150   = 5.0599e-6;
    RDouble rho150 = 2.0752e-9;
    RDouble t160   = 1022.20;
    RDouble p160   = 3.6929e-6;
    RDouble rho160 = 1.2336e-9;
    RDouble t170   = 1103.40;
    RDouble p170   = 2.7915e-6;
    RDouble rho170 = 7.8155e-10;
    RDouble t190   = 1205.40;
    RDouble p190   = 1.6845e-6;
    RDouble rho190 = 3.5807e-10;
    RDouble t230   = 1322.30;
    RDouble p230   = 6.7138e-7;
    RDouble rho230 = 1.5640e-10;
    RDouble t300   = 1432.10;
    RDouble p300   = 1.8828e-7;
    RDouble rho300 = 1.9159e-11;
    RDouble t400   = 1487.40;
    RDouble p400   = 4.0278e-8;
    RDouble rho400 = 2.8028e-12;
    RDouble t500   = 1499.20;
    RDouble p500   = 1.0949e-8;
    RDouble rho500 = 5.2148e-13;
    RDouble t600   = 1506.10;
    RDouble p600   = 3.4475e-9;
    RDouble rho600 = 1.1367e-13;
    RDouble t700   = 1507.60;
    RDouble p700   = 1.1908e-9;
    RDouble rho700 = 1.5270e-13;

    RDouble rho = 0.0;

    if (h <= 11019.0) 
    {
        RDouble al1 = ( t11 - t0 )/11019.0;
        t   = t0 + al1 * h;
        p   = p0 * pow( t/t0, -g/(r * al1 ) );
        rho = rho0 * pow( t/t0, -1.0 - g/(r*al1) );
    }
    else if (h <= 20063.0)
    {
        t   = t11;
        p   = p11   * exp( -g*(h-11019.0)/(r*t11) );
        rho = rho11 * exp( -g*(h-11019.0)/(r*t11) );
    }
    else if (h <= 32162.0)
    {
        RDouble al2 = (t32-t20)/(32162.0-20063.0);
        t   = t11+al2*(h-20063.0);
        p   = p20 * pow(t/t11, -g/(r*al2));
        rho = rho20 * pow(t/t11, -1.0-g/(r*al2));
    }
    else if (h <= 47350.0)
    {
        RDouble al3 = (t47-t32)/(47350.0-32162.0);
        t   = t32+al3*(h-32162.0);
        p   = p32 * pow(t/t32, -g/(r*al3));
        rho = rho32 * pow(t/t32, -1.0-g/(r*al3));
    }

    else if (h <= 52429.0)
    {
        t   = t47;
        p   = p47*exp(-g*(h-47350.0)/(r*t47));
        rho = rho47*exp(-g*(h-47350.0)/(r*t47));
    }
    else if(h <= 61591.0)
    {
        RDouble al4 =( t61-t52)/(61591.0-52429.0);
        t = t47+al4*(h-52429.0);
        p = p52 * pow(t/t47, -g/(r*al4));
        rho = rho52* pow(t/t47,-1.0-g/(r*al4));
    }
    else if(h <= 79994.0)
    {
        RDouble al5=(t79-t61)/(79994.0-61591.0);
        t=t61+al5*(h-61591.0);
        p=p61*pow(t/t61, -g/(r*al5));
        rho=rho61*pow(t/t61,-1.0-g/(r*al5));
    }
    else if(h <= 90000.0)
    {
        t=t79;
        p=p79*exp(-g*(h-79994.0)/(r*t79));
        rho=rho79*exp(-g*(h-79994.0)/(r*t79));
    }
    else if(h <= 100000.0)
    {
        RDouble  al6=(t100-t90)/10000.0;
        t=t79+al6*(h-90000.0);
        p=p90* pow(t/t79,-g/(r*al6));
        rho=rho90* pow(t/t79,-1.0-g/(r*al6));
    }
    else if(h <= 110000.0)
    {
        RDouble al7=(t110-t100)/10000.0;
        t=t100+al7*(h-100000.0);
        p=p100*pow(t/t100, -g/(r*al7));
        rho=rho100*(t/t100, -1.0-g/(r*al7));
    }
    else if(h <= 120000.0)
    {
        RDouble al8=(t120-t110)/10000.0;
        t=t110+al8*(h-110000.0);
        p=p110*pow(t/t110, -g/(r*al8));
        rho=rho110*pow(t/t110, -1.0-g/(r*al8));
    }
    else if(h <= 150000.0)
    {
        RDouble al9=(t150-t120)/30000.0;
        t=t120+al9*(h-120000.0);
        p=p120*pow(t/t120,-g/(r*al9));
        rho=rho120*pow(t/t120,-1.0-g/(r*al9));
    }
    else if(h <= 160000.0)
    {
        RDouble al10=(t160-t150)/10000.0;
        t=t150+al10*(h-150000.0);
        p=p150*pow(t/t150, -g/(r*al10));
        rho=rho150*pow(t/t150, -1.0-g/(r*al10));
    }
    else if(h <= 170000.0)
    {
        RDouble al11=(t170-t160)/10000.0;
        t=t160+al11*(h-160000.0);
        p=p160*pow(t/t160, -g/(r*al11));
        rho=rho160*pow(t/t160,-1.0-g/(r*al11));
    }
    else if(h <= 190000.0)
    {
        RDouble al12=(t190-t170)/20000.0;
        t=t170+al12*(h-170000.0);
        p=p170*pow(t/t170, -g/(r*al12));
        rho=rho170*pow(t/t170,-1.0-g/(r*al12));
    }
    else if(h <= 230000.0)
    {
        RDouble al13=(t230-t190)/40000.0;
        t=t190+al13*(h-190000.0);
        p=p190*pow(t/t190,-g/(r*al13));
        rho=rho190*pow(t/t190,-1.0-g/(r*al13));
    }
    else if(h <= 300000.0)
    {
        RDouble al14=(t300-t230)/70000.0;
        t=t230+al14*(h-230000.0);
        p=p230*pow(t/t230, -g/(r*al14));
        rho=rho230*pow(t/t230, -1.0-g/(r*al14));
    }
    else if(h <= 400000.0)
    {
        RDouble al15=(t400-t300)/100000.0;
        t=t300+al15*(h-300000.0);
        p=p300*pow(t/t300, -g/(r*al15));
        rho=rho300*pow(t/t300, -1.0-g/(r*al15));
    }
    else if(h <= 500000.0)
    {
        RDouble al16=(t500-t400)/100000.0;
        t=t400+al16*(h-400000.0);
        p=p400*pow(t/t400, -g/(r*al16));
        rho=rho400*pow(t/t400, -1.0-g/(r*al16));
    }
    else if(h <= 600000.0)
    {
        RDouble al17=(t600-t500)/100000.0;
        t=t500+al17*(h-500000.0);
        p=p500*pow(t/t500, -g/(r*al17));
        rho=rho500*pow(t/t500, -1.0-g/(r*al17));
    }
    else if(h <= 700000.0)
    {
        RDouble al18=(t700-t600)/100000.0;
        t=t600+al18*(h-600000.0);
        p=p600*pow(t/t600, -g/(r*al18));
        rho=rho600*pow(t/t600, -1.0-g/(r*al18));
    }
    else
    {
        TK_Exit::ExceptionExit("Error: height must be less than 700km!!!\n", true);
    }

    a   = sqrt(1.4*r*t);
    p   = p * 100.0;
    den = rho;
}

void MixedGas::ComputeCoefficientInMXDQ(RDouble *primitiveVars, RDouble &alpha, RDouble &beta, RDouble *speciesBeta)
{
    //! The arguments alpha, beta and beta_s are referred to the formular (A.5) in the appendix A of the PHengLEI Theory manual.
    //! alpha = R/(M*cv)
    //! beta = R/(M*cv)*(0.5*V^2 - ens) + R*T/Mns
    //! betas =  R/(M*cv) * (ens - es) + R*T(1/Ms - 1/Mns)
    int nSpecies = this->GetNumberOfSpecies();
    for (int s = 0; s < nSpecies; ++ s)
    {
        massFractions[s] = primitiveVars[this->nm + s];
    }

    //! Obtain the primitive values.
    using namespace IDX;
    RDouble &density =   primitiveVars[IR];
    RDouble &uVelocity = primitiveVars[IU];
    RDouble &vVelocity = primitiveVars[IV];
    RDouble &wVelocity = primitiveVars[IW];
    RDouble &pressure =  primitiveVars[IP];
    RDouble squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;
    //! Obtain the non-dimensional universal gas constant.
    RDouble gasConstant = this->GetUniversalGasConstant();
    //! Obtain the molecular weight of mixture gas.
    RDouble mixedGasMass = this->GetMixedGasMolecularWeight(massFractions);
    //! Compute the temperature using the state equation of mixture gas.
    RDouble temperature = pressure * mixedGasMass / (density * gasConstant);

    //! Compute the internal energy of last species with the label ns-1.
    this->ComputeSpeciesConstantVolumeSpecificHeat(temperature, speciesCv);
    //this->ComputeSpeciesEnthalpy(temperature, speciesCv, speciesEnergyTr, speciesEnthalpy);
    RDouble mixedGasCv = 0.0;
    this->ComputeMixtureByMassFraction(massFractions, speciesCv, mixedGasCv);
    //! Compute the first coefficient called alpha.
    alpha = gasConstant / (mixedGasMass * mixedGasCv);

    //! Obtain the molecular weight of each species.
    RDouble *speciesMassReciprocal = this->molecularProperty->GetMolecularWeightReciprocal();
    RDouble nitrogenEnergy = speciesCv[nLastSpeciesIndex] * temperature;//speciesEnergyTr[nLastSpeciesIndex];
    RDouble electronEnergy = 0.0, electronMassReciprocal = 1.0 / electronMass;
    if (nElectronIndex >= 0)
    {
        electronEnergy = speciesCv[nElectronIndex] * temperature;
        //electronEnergy = speciesEnergyTr[nElectronIndex];
        electronMassReciprocal = speciesMassReciprocal[nElectronIndex];
    }

    //! Compute the second coefficient called beta.
    RDouble alpha1 = 0.5 * squareVelocity - nitrogenEnergy;
    RDouble alpha2 = speciesMassReciprocal[nLastSpeciesIndex];
    beta = alpha * alpha1 + gasConstant * temperature * alpha2;

    //! Compute the other coefficient stored in speciesBeta.
    RDouble speciesGama, speciesEta, tempCoef;
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    speciesBeta[nLastSpeciesIndex] = 0.0;
    for (int s = 0; s < nLastSpeciesIndex; ++ s)
    {
        tempCoef = ionType[s] * speciesMassReciprocal[s] / electronMassReciprocal;

        speciesGama = speciesMassReciprocal[s] - speciesMassReciprocal[nLastSpeciesIndex];
        speciesGama += tempCoef * (electronMassReciprocal - speciesMassReciprocal[nLastSpeciesIndex]);

        speciesEta = nitrogenEnergy - speciesCv[s] * temperature;
        //speciesEta = nitrogenEnergy - speciesEnergyTr[s];
        speciesEta += tempCoef * (nitrogenEnergy - electronEnergy);

        speciesBeta[s] = alpha * speciesEta + gasConstant * temperature * speciesGama;
    }
}

void MixedGas::ComputeCoefficientInMXDQR(RDouble *primitiveVars, RDouble trTemperature, RDouble squareVelocity, RDouble *speciesCvs, RDouble totalCv, RDouble &alpha, RDouble &beta, RDouble *speciesBeta)
{
    //! The arguments alpha, beta and beta_s are referred to the formular (A.5) in the appendix A of the PHengLEI Theory manual.
    //! alpha = R/(M*cv)
    //! beta = R/(M*cv)*(0.5*V^2 - ens) + R*T/Mns
    //! betas =  R/(M*cv) * (ens - es) + R*T(1/Ms - 1/Mns)
    //! Obtain the primitive values.

    RDouble gasConstant = this->GetUniversalGasConstant();
    RDouble *speciesMassReciprocal = this->molecularProperty->GetMolecularWeightReciprocal();
    int *ionType = molecularProperty->GetIonTypeOfSpecies();

    RDouble mixedGasMass = 0.0;
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVars[this->nm + i];
        mixedGasMass += massFractions[i] * speciesMassReciprocal[i];
    }

    alpha = gasConstant * mixedGasMass / totalCv;

    //! Obtain the molecular weight of each species.

    RDouble nitrogenEnergy = speciesCvs[nLastSpeciesIndex] * trTemperature;   //speciesEnergyTr[nLastSpeciesIndex];

    RDouble electronEnergy = 0.0, electronMassReciprocal = 1.0 / electronMass;
    if (nElectronIndex >= 0)
    {
        electronEnergy = speciesCvs[nElectronIndex] * trTemperature;
        //electronEnergy = speciesEnergyTr[nElectronIndex];
        electronMassReciprocal = speciesMassReciprocal[nElectronIndex];
    }

    //! Compute the second coefficient called beta.
    RDouble alpha1 = 0.5 * squareVelocity - nitrogenEnergy;
    RDouble alpha2 = speciesMassReciprocal[nLastSpeciesIndex];
    beta = alpha * alpha1 + gasConstant * trTemperature * alpha2;

    //! Compute the other coefficient stored in speciesBeta.
    RDouble speciesGama, speciesEta, tempCoef;

    for (int s = 0; s < nLastSpeciesIndex; ++ s)
    {
        tempCoef = ionType[s] * speciesMassReciprocal[s] / electronMassReciprocal;

        speciesGama = speciesMassReciprocal[s] - speciesMassReciprocal[nLastSpeciesIndex];
        speciesGama += tempCoef * (electronMassReciprocal - speciesMassReciprocal[nLastSpeciesIndex]);

        speciesEta = nitrogenEnergy - speciesCvs[s] * trTemperature;
        //speciesEta = nitrogenEnergy - speciesEnergyTr[s];
        speciesEta += tempCoef * (nitrogenEnergy - electronEnergy);

        speciesBeta[s] = alpha * speciesEta + gasConstant * trTemperature * speciesGama;
    }
    speciesBeta[nLastSpeciesIndex] = 0.0;
}

//! To obtain the coefficients of computing the M*dQ in the function MVXDQ_STD_LP(), added by LiPeng on Jan. 21, 2019.
void MixedGas::ComputeCoefficientInMVXDQ(RDouble *prim, RDouble K, RDouble rD, RDouble *theta_s, RDouble *phi_s, RDouble *Ds)
{
    //! K=visl/Prl + vist/Prt denotes the thermal conductivity coefficient, rD=visl/Scl + vist/Sct is the species diffusion coefficient,\n
    //! Ds is an array of saving the diffusion coefficient of each species.

    int ns = this->GetNumberOfSpecies();
    //! To obtain the mass fraction of each species.
    RDouble *fs = new RDouble[ns];
    RDouble *xs = new RDouble[ns];
    for (int i = 0; i < ns; i++)
    {
        fs[i] = prim[i + this->nm];
    }

    //! Obtain the mole fraction of each species.
    this->ComputeMoleFractionByMassFraction(fs, xs);

    //! Obtain the density and pressure.
    RDouble &ro = prim[0];
    RDouble &pm = prim[4];
    RDouble R = this->GetUniversalGasConstant();    //! The non-dimensional universal gas constant.
    RDouble M = this->GetMixedGasMolecularWeight(fs);     //! The molecular weight of mixture gas.
    RDouble tm = pm * M / (ro * R);     //! Compute the temperature using the state equation of mixture gas.

    RDouble cp;    //! Save the value of constant pressure specific heat.
    this->ComputeConstantPressureSpecificHeatByMassFraction(fs,tm,cp);
    RDouble lambda = cp * K;      //! K = visl/Prl + vist/Prt denotes the thermal conductivity coefficient.

    //! Compute the diffusion coefficient of each species.
    for (int i = 0; i < ns; i++)
    {
        Ds[i] = (rD/ro) * (1.0 - fs[i])/(1.0 - xs[i]);     //! rD=visl/Scl + vist/Sct.
    }

    RDouble cv;//!Save the value of constant volume specific heat.
    this->ComputeConstantVolumeSpecificHeatByMassFraction(fs, tm, cv);
    RDouble e = cv * tm;//!the total internal energy of mixture gas

    //! Compute the internal energy of last species with label ns-1 in the collection.
    RDouble *cvs = new RDouble[ns];    //! Save the value of constant volume specific heat of each species.
    this->ComputeSpeciesConstantVolumeSpecificHeat(tm,cvs);
    RDouble ens = cvs[ns-1] * tm;
    RDouble *Ms = this->molecularProperty->GetMolecularWeight();    //! Obtain the molecular weight of each species.

    //! Compute the coefficient called theta_s.
    RDouble ei = 0.0;
    for (int i = 0; i < ns; i++)
    {
        ei = cvs[i] * tm;
        theta_s[i] = (ei - ens) - M * (1.0/Ms[i] - 1.0/Ms[ns-1]) * e;
    }

    //! Compute the coefficient called phi_s.
    RDouble *cps = new RDouble[ns];    //! Save the value of constant pressure specific heat of each species.
    this->ComputeSpeciesConstantPressureSpecificHeat(tm, cps);
    RDouble roDns = ro * Ds[ns-1] * cps[ns-1] * tm;
    RDouble hi = 0.0;

    //! Compute the coefficient called phi_s of each species.
    for (int i = 0; i < ns; i++)
    {//!lambda=cp*(visl/Prl + vist/Prt)
        hi = cps[i] * tm;
        phi_s[i] = -lambda * M * (1.0/Ms[i] - 1.0/Ms[ns-1]) * tm + (ro * Ds[i] * hi  - roDns);
    }

    //! Remove the dynamic memory space.
    delete [] fs;
    delete [] xs;
    delete [] cvs;
    delete [] cps;
}

void MixedGas::GetEarthFullyCatalyticMassFraction(RDouble *fs, RDouble *fcw)
{
    RDouble *Ms = this->molecularProperty->GetMolecularWeightDimensional();
    RDouble M = 0.0;
    //! Compute the total molecular weight.
    for (int i = 0; i < this->numberOfSpecies; i++)
    {
        fcw[i] = 0.0; //Initialization.
        M += fs[i] / Ms[i];
    }
    M = 1.0 / M;

    RDouble *xs = new RDouble[this->numberOfSpecies];
    //! Compute mole fraction of each species.
    for (int i = 0; i < this->numberOfSpecies; i++)
    {
        xs[i] = fs[i] * M / Ms[i];
    }

    //! To compute the mass of oxygen.
    RDouble co2 = 0.0, cn2 = 0.0;
    int coef = 0;
    RDouble Mo = 15.9994E-3;     //! The molecular weight.
    string *species_name = this->molecularProperty->GetNameOfSpecies();
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        if (FindChemicalElement(species_name[i], 'O', coef))
        {
            co2 += xs[i] * coef * Mo;
        }
    }
    delete []xs;
    //Compute the mass fraction of Oxygen.
    co2 = co2 / M;
    //Compute the mass fraction of Nitrogen.
    cn2 = MAX(1.0 - co2, 0.0);

    int nIndex = this->GetSpeciesIndex("O2");
    if (nIndex >= 0)
    {
        fcw[nIndex] = co2;
    }
    nIndex = this->GetSpeciesIndex("N2");
    if (nIndex >= 0)
    {
        fcw[nIndex] = cn2;
    }
}

void MixedGas::GetMarsFullyCatalyticMassFraction(RDouble *fs, RDouble *fcw)
{
    RDouble *Ms = this->molecularProperty->GetMolecularWeightDimensional();
    RDouble M = 0.0;
    //! Compute the total molecular weight.
    for (int i = 0; i < this->numberOfSpecies; i++)
    {
        fcw[i] = 0.0; //Initialization.
        M += fs[i] / Ms[i];
    }
    M = 1.0 / M;

    RDouble *xs = new RDouble[this->numberOfSpecies];
    //! Compute mole fraction of each species.
    for (int i = 0; i < this->numberOfSpecies; i++)
    {
        xs[i] = fs[i] * M / Ms[i];
    }

    //! To compute the mass of Oxygen, Nitrogen and Carbon dioxide.
    RDouble cO2 = 0.0, cN2 = 0.0, cCO2 = 0.0;
    int coef = 0;
    RDouble atomicNitrogenMass = 1.400674E-2, carbonDioxideMass = 4.401E-2; //! The molecular weight(kg/mol).
    string *species_name = this->molecularProperty->GetNameOfSpecies();

    RDouble moleFracationCO2 = 0.0;
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        //Compute the total mass of Nitrogen.
        if (FindChemicalElement(species_name[i], 'N', coef))
        {
            cN2 += xs[i] * coef * atomicNitrogenMass;
        }
        //Compute the total mass of Carbon.
        if (FindChemicalElement(species_name[i], 'C', coef))
        {
            moleFracationCO2 += xs[i] * coef; //Obtain the total mole fractions of CO2;
        }
    }
    delete []xs;

    //Compute the mass fraction of Nitrogen.
    cN2 = cN2 / M;
    //Compute the mass fraction of Carbon dioxide.
    cCO2 = moleFracationCO2 * carbonDioxideMass / M;
    //Compute the mass fraction of Oxygen.
    cO2 = MAX(1.0 - cN2 - cCO2, 0.0);

    if(cCO2 > 1.0 && cN2 == 0 && cO2 == 0)
    {
        cCO2 = 1.0;
    }

    int nIndex = this->GetSpeciesIndex("O2");
    if (nIndex >= 0)
    {
        fcw[nIndex] = cO2;
    }
    nIndex = this->GetSpeciesIndex("N2");
    if (nIndex >= 0)
    {
        fcw[nIndex] = cN2;
    }
    nIndex = this->GetSpeciesIndex("CO2");
    if (nIndex >= 0)
    {
        fcw[nIndex] = cCO2;
    }
}

bool MixedGas::FindChemicalElement(string species_name, char element_name, int &coef)
{
    bool bflag = false;
    coef = 0;
    int i = static_cast<int> (species_name.find(element_name));
    if (i >= 0)
    {
        coef = 1;
        bflag = true;
        int n = static_cast<int> (species_name.length());
        if (i + 1 < n)
        {
            if (species_name[i + 1] >= '1' && species_name[i + 1] <= '9')
            {
                coef = atoi(species_name.substr(i + 1, 1).c_str());
            }

            if (element_name == 'N' && species_name[i + 1] == 'a')  //Na
            {
                coef = 0.0;
                bflag = false;   //Õâ¸öµØ·½Ã»¿¼ÂÇNaµÄµª»¯Îï
            }

            if (element_name == 'S' && species_name[i + 1] == 'i')  //Si
            {
                coef = 0.0;     //
                bflag = false;   //Õâ¸öµØ·½Ã»¿¼ÂÇSiµÄÁò»¯Îï
            }
        }
    }

    return bflag;
}

int MixedGas::GetSpeciesIndex(string species_name)
{
    string *name_array = this->molecularProperty->GetNameOfSpecies();
    for (int i = 0; i < this->numberOfSpecies; i++)
    {
        if (name_array[i] == species_name)
        {
            return i;
        }
    }
    return -1;
}

void MixedGas::ComputeTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *speciesTransRotationCv)
{
    //! The dimensional value of translation-rotation specific heat at constant volume.
    ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesTransRotationCv);

    //! Translate to the non-dimensional values.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesTransRotationCv[i] *= referenceSpecificHeat;
    }
}

void MixedGas::ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *speciesTransRotationCv)
{
    RDouble gasConstant = rjmk; //Unit: J/(kg*K)
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();
    RDouble coef = 1.0;
    for (int i = 0; i <= this->nLastSpeciesIndex; ++ i)
    {
        speciesTransRotationCv[i] = 1.5 * gasConstant / speciesMass[i];
        coef = (RDouble)(thermodynamicProperties[i].nType - 1);
        //! Polyatom species.
        speciesTransRotationCv[i] += coef * gasConstant / speciesMass[i];
    }
    if (this->nElectronIndex >= 0)
    {
        speciesTransRotationCv[this->nElectronIndex] = 0.0;
    }
}

void MixedGas::ComputeTranslationAndRotationEnergy(RDouble transRotationTemperature, RDouble *speciesTransRotationEnergy)
{
    //! Obtain the translation-rotation specific heat at constant volume.
    ComputeTranslationAndRotationSpecificHeatAtConstantVolume(speciesTransRotationEnergy);

    //! Obtain the thermodynamic properties of species.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();

    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesTransRotationEnergy[i]  = speciesTransRotationEnergy[i] * transRotationTemperature;
        speciesTransRotationEnergy[i] += thermodynamicProperties[i].hs0 / referenceEnergy;
    }
}

void MixedGas::ComputeVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *speciesVibrationCv)
{
    //! The dimensional value of vibration specific heat at constant volume.
    ComputeDimensionalVibrationSpecificHeatAtConstantVolume(vibrationTemperature, speciesVibrationCv);

    //! Get non-dimensional values.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesVibrationCv[i] *= referenceSpecificHeat;
    }
}

void MixedGas::ComputeDimensionalVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *speciesVibrationCv)
{
    RDouble gasConstant = rjmk; //Unit: J/(kg*K)
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    //! The reference temperature and the dimensional vibration temperature.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble Tv = vibrationTemperature * referenceTemperature;

    if (nTEnergyModel >= 1)
    {
        RDouble speciesCvtr[MAX_SPECIES_NUM], speciesCvv[MAX_SPECIES_NUM];
        if (nTEnergyModel == 1)
        {
            this->ComputeSpecificHeatByFitMethod(vibrationTemperature, speciesCvv, true);
        }
        else
        {
            this->ComputeSpecificHeatByPolynomialFitMethod(vibrationTemperature, speciesCvv, true);
        }
        this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCvtr);
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType < 2)
            {
                speciesVibrationCv[i] = 0.0;
                continue;
            }
            //! Polyatom species.
            speciesVibrationCv[i] = speciesCvv[i] - speciesCvtr[i];
        }
    }
    else
    {
        RDouble theta = 0.0, expValue = 0.0, powValue = 0.0;
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType < 2)
            {
                speciesVibrationCv[i] = 0.0;
                continue;
            }
            //! Polyatom species.
            speciesVibrationCv[i] = 0.0;
            for (int m = 0; m < thermodynamicProperties[i].Tvs.nMode; ++ m)
            {
                theta = thermodynamicProperties[i].Tvs.Tv[m] / (Tv + SMALL);

                expValue = MIN(exp(theta), LARGE);
                powValue = theta / ((expValue - 1.0) + SMALL);
                powValue = powValue * powValue;

                speciesVibrationCv[i] += thermodynamicProperties[i].Tvs.g[m] * (gasConstant / speciesMass[i]) * expValue * powValue;
            }
        }
    }
}

void MixedGas::ComputeDimensionalVibrationEnergy(RDouble vibrationTemperature, RDouble *speciesVibrationEnergy)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();

    //! Obtain the reference temperature and velocity.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble Tv = vibrationTemperature * referenceTemperature;

    if (nTEnergyModel >= 1)
    {
        RDouble speciesCvtr[MAX_SPECIES_NUM], speciesEnthaly[MAX_SPECIES_NUM];
        if (nTEnergyModel == 1)
        {
            this->ComputeEnthalpyByFitMethod(vibrationTemperature, speciesEnthaly, true);
        }
        else
        {
            this->ComputeEnthalpyByPolynomialFitMethod(vibrationTemperature, speciesEnthaly, true);
        }
        this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCvtr);
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType < 2)
            {
                speciesVibrationEnergy[i] = 0.0;
                continue;
            }
            //! Polyatom species.
            speciesVibrationEnergy[i] = speciesEnthaly[i] - gasConstant * Tv / speciesMass[i];
            speciesVibrationEnergy[i] -= speciesCvtr[i] * Tv + thermodynamicProperties[i].hs0;
        }
    }
    else
    {
        RDouble theta = 0.0, expValue = 0.0, powValue = 0.0;
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType < 2)
            {
                speciesVibrationEnergy[i] = 0.0;
                continue;
            }
            //! Polyatom species.
            speciesVibrationEnergy[i] = 0.0;
            for (int m = 0; m < thermodynamicProperties[i].Tvs.nMode; ++ m)
            {
                theta = thermodynamicProperties[i].Tvs.Tv[m] / (Tv + SMALL);
                expValue = MIN(exp(theta), LARGE);
                powValue = thermodynamicProperties[i].Tvs.Tv[m] / ((expValue - 1.0) + SMALL);

                speciesVibrationEnergy[i] += thermodynamicProperties[i].Tvs.g[m] * (gasConstant / speciesMass[i]) * powValue;
            }
        }
    }
}

void MixedGas::ComputeVibrationEnergy(RDouble vibrationTemperature, RDouble *speciesVibrationEnergy)
{
    //! Obtain the dimensional value of vibration energy.
    ComputeDimensionalVibrationEnergy(vibrationTemperature, speciesVibrationEnergy);

    //! Translate to non-dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesVibrationEnergy[i] /= referenceEnergy;
    }
}

void MixedGas::ComputeElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *speciesElectronCv)
{
    //! The dimensional value of electron specific heat at constant volume.
    ComputeDimensionalElectronSpecificHeatAtConstantVolume(electronTemperature, speciesElectronCv);

    //! Translate to the non-dimensional values.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesElectronCv[i] *= referenceSpecificHeat;
    }
}

void MixedGas::ComputeDimensionalElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *speciesElectronCv)
{
    RDouble gasConstant = rjmk;      //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();

    //! Obtain the reference temperature and compute the dimensional value of electron temperature.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble Te = electronTemperature * referenceTemperature;

    if (nTEnergyModel >= 1)
    {
        RDouble speciesCvtr[MAX_SPECIES_NUM], speciesCvv[MAX_SPECIES_NUM];
        if (nTEnergyModel == 1)
        {
            this->ComputeSpecificHeatByFitMethod(electronTemperature, speciesCvv, false);
        }
        else
        {
            this->ComputeSpecificHeatByPolynomialFitMethod(electronTemperature, speciesCvv, false);
        }
        this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCvtr);
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType >= 2)
            {
                speciesElectronCv[i] = 0.0;
                continue;
            }
            speciesElectronCv[i] = speciesCvv[i] - speciesCvtr[i];
        }
    }
    else
    {
        RDouble theta = 0.0, expValue = 0.0, gs0 = 0.0, gs1 = 0.0, powValue = 0.0;
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType == 0) //electron
            {
                speciesElectronCv[i] = 1.5 * gasConstant / speciesMass[i];
                continue;
            }

            //! The degeneracy of species.
            gs0 = thermodynamicProperties[i].gs0;
            gs1 = thermodynamicProperties[i].gs1;

            //! The monoaton species and polyatom species.
            theta = thermodynamicProperties[i].Tes / (Te + SMALL);
            expValue = MIN(exp(theta), LARGE);
            powValue = theta / (gs0 * expValue + gs1);
            powValue = powValue * powValue;

            speciesElectronCv[i] = (gasConstant / speciesMass[i]) * (gs0 * gs1) * expValue * powValue;
        }
    }
}

void MixedGas::ComputeDimensionalElectronEnergy(RDouble electronTemperature, RDouble *speciesElectronEnergy)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();

    //! Obtain the reference temperature and compute the dimensional value of electron temperature.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble Te = electronTemperature * referenceTemperature;

    if (nTEnergyModel >= 1)
    {
        RDouble speciesCvtr[MAX_SPECIES_NUM], speciesEnthaly[MAX_SPECIES_NUM];
        if (nTEnergyModel == 1)
        {
            this->ComputeEnthalpyByFitMethod(electronTemperature, speciesEnthaly, false);
        }
        else
        {
            this->ComputeEnthalpyByPolynomialFitMethod(electronTemperature, speciesEnthaly, false);
        }
        this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCvtr);
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType >= 2)
            {
                speciesElectronEnergy[i] = 0.0;
                continue;
            }
            speciesElectronEnergy[i] = speciesEnthaly[i] - gasConstant * Te / speciesMass[i];
            speciesElectronEnergy[i] -= speciesCvtr[i] * Te + thermodynamicProperties[i].hs0;
        }
    }
    else
    {
        RDouble theta = 0.0, expValue = 0.0, gs0 = 0.0, gs1 = 0.0;
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType == 0) //electron
            {
                speciesElectronEnergy[i] = 1.5 * gasConstant * Te / speciesMass[i];
                continue;
            }

            //! The monoaton species and polyatom species.
            theta = thermodynamicProperties[i].Tes / (Te + SMALL);
            //! The degeneracy of species.
            gs0 = thermodynamicProperties[i].gs0;
            gs1 = thermodynamicProperties[i].gs1;
            expValue = MIN(exp(theta), LARGE);

            speciesElectronEnergy[i] = (gasConstant / speciesMass[i]) * thermodynamicProperties[i].Tes * (gs1 / (gs0 * expValue + gs1));
        }
    }
}

void MixedGas::ComputeElectronEnergy(RDouble electronTemperature, RDouble *speciesElectronEnergy)
{
    //! Obtain the dimensional value of electron energy.
    ComputeDimensionalElectronEnergy(electronTemperature, speciesElectronEnergy);

    //! Translate to non-dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesElectronEnergy[i] /= referenceEnergy;
    }
}

void MixedGas::ComputeSpeciesEnthalpy(RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble *speciesEnthalpy)
{
    //! Compute energies.
    ComputeTranslationAndRotationEnergy(transRotationTemperature, this->speciesEnergyTr);
    ComputeVibrationEnergy(vibrationTemperature, this->speciesEnergyVb);
    ComputeElectronEnergy(electronTemperature, this->speciesEnergyEe);

    RDouble gasConstant = GetUniversalGasConstant();
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();
    if (this->nElectronIndex >= 0)
    {
        //! Electron.
        speciesEnthalpy[nElectronIndex] = 2.5 * gasConstant * electronTemperature / speciesMass[nElectronIndex];
    }
    for (int i = 0; i <= nLastSpeciesIndex; ++ i)
    {
        //! The monoaton species and polyatom species.
        speciesEnthalpy[i] = speciesEnergyTr[i] + speciesEnergyVb[i] + speciesEnergyEe[i];
        speciesEnthalpy[i] += gasConstant * transRotationTemperature / speciesMass[i];
    }
}

RDouble MixedGas::GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *massFraction)
{
    //! Get dimensional value of translation-rotation specific heat at constant volume.
    RDouble transRotationCv = GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(massFraction);

    //! Translate to the non-dimensional values.
    return (transRotationCv * referenceSpecificHeat);
}

RDouble MixedGas::GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *massFraction)
{
    ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);

    RDouble transRotationCv = 0.0;     //! The translation-rotation specific heat at constant volume.
    ComputeMixtureByMassFraction(massFraction, speciesCv, transRotationCv);

    return transRotationCv;
}

RDouble MixedGas::GetMixedGasVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *massFraction)
{
    //! Get dimensional value of vibration specific heat at constant volume.
    RDouble vibrationCv = GetMixedGasDimensionalVibrationSpecificHeatAtConstantVolume(vibrationTemperature, massFraction);

    // Translate to the non-dimensional values.
    return (vibrationCv * referenceSpecificHeat);
}

RDouble MixedGas::GetMixedGasDimensionalVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *massFraction)
{
    ComputeDimensionalVibrationSpecificHeatAtConstantVolume(vibrationTemperature, speciesCv);

    //! The vibration specific heat at constant volume.
    RDouble vibrationCv = 0.0;
    ComputeMixtureByMassFraction(massFraction, speciesCv, vibrationCv);

    return vibrationCv;
}

RDouble MixedGas::GetMixedGasElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *massFraction)
{
    //! Get dimensional value of electron specific heat at constant volume.
    RDouble electronCv = GetMixedGasDimensionalElectronSpecificHeatAtConstantVolume(electronTemperature, massFraction);

    //! Translate to the non-dimensional values.
    return (electronCv * referenceSpecificHeat);
}

RDouble MixedGas::GetMixedGasDimensionalElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *massFraction)
{
    this->ComputeDimensionalElectronSpecificHeatAtConstantVolume(electronTemperature, speciesCv);

    RDouble electronCv = 0.0;     //! The electron specific heat at constant volume.
    ComputeMixtureByMassFraction(massFraction, speciesCv, electronCv);

    return electronCv;
}

RDouble MixedGas::GetMixedGasTranslationAndRotationEnergy(RDouble transRotationTemperature, RDouble *massFraction)
{
    //! To obtain the translation-rotation energy.
    ComputeTranslationAndRotationEnergy(transRotationTemperature, speciesEnergyTr);

    RDouble transRotationEnergy = 0.0;
    ComputeMixtureByMassFraction(massFraction, speciesEnergyTr, transRotationEnergy);

    return transRotationEnergy;
}

RDouble MixedGas::GetMixedGasDimensionalVibrationEnergy(RDouble vibrationTemperature, RDouble *massFraction)
{
    //! To obtain the vibration energy.
    this->ComputeDimensionalVibrationEnergy(vibrationTemperature, speciesEnergyVb);

    RDouble vibrationEnergy = 0.0;
    ComputeMixtureByMassFraction(massFraction, speciesEnergyVb, vibrationEnergy);

    return vibrationEnergy;
}

RDouble MixedGas::GetMixedGasVibrationEnergy(RDouble vibrationTemperature, RDouble *massFraction)
{
    RDouble vibrationEnergy = GetMixedGasDimensionalVibrationEnergy(vibrationTemperature, massFraction);
    //! Translate to non-dimensional value.
    vibrationEnergy = vibrationEnergy / referenceEnergy;

    return vibrationEnergy;
}

RDouble MixedGas::GetMixedGasDimensionalElectronEnergy(RDouble electronTemperature, RDouble *massFraction)
{
    //! To obtain the electron energy.
    this->ComputeDimensionalElectronEnergy(electronTemperature, speciesEnergyEe);

    RDouble electronEnergy = 0.0;
    ComputeMixtureByMassFraction(massFraction, speciesEnergyEe, electronEnergy);

    return electronEnergy;
}

RDouble MixedGas::GetMixedGasElectronEnergy(RDouble electronTemperature, RDouble *massFraction)
{
    RDouble electronEnergy = GetMixedGasDimensionalElectronEnergy(electronTemperature, massFraction);
    //! Translate to non-dimensional value.
    electronEnergy = electronEnergy / referenceEnergy;

    return electronEnergy;
}

RDouble MixedGas::ComputeMixedGasEnthalpy(RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble *massFraction)
{
    //! To obtain the enthalpies of species.
    ComputeSpeciesEnthalpy(transRotationTemperature, vibrationTemperature, electronTemperature, speciesEnthalpy);

    RDouble enthalpy = 0.0;
    ComputeMixtureByMassFraction(massFraction, speciesEnthalpy, enthalpy);

    return enthalpy;
}

RDouble MixedGas::ComputeMixedGasEnthalpy(RDouble *primitiveVars, RDouble *temperatures)
{
    using namespace IDX;
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        massFractions[s] = primitiveVars[this->nm + s];
    }
    RDouble staticEnthalpy = 0.0;
    if (this->ntmodel == 1)
    {
        staticEnthalpy = this->GetMixtureGasEnthalpy(massFractions, temperatures[ITT]);
     }
    else
    {
        staticEnthalpy = this->ComputeMixedGasEnthalpy(temperatures[ITT], temperatures[ITV], temperatures[ITE], massFractions);
    }
    return staticEnthalpy;
}

RDouble MixedGas::GetTotalFormationEnthalpy(RDouble *massFraction)
{
    //! To obtain the thermodynamic properties of species.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();

    //! To record the total enthalpy of formation.
    RDouble formationEnthalpy = 0.0;
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        formationEnthalpy += massFraction[i] * thermodynamicProperties[i].hs0;
    }
    return (formationEnthalpy / referenceEnergy);
}

RDouble MixedGas::ComputeTranslationAndRotationTemperatureByPrimitiveVariables(RDouble *primitiveVariables, RDouble electronPressure/* = 0.0*/)
{
    using namespace IDX;
    RDouble density = primitiveVariables[IR];
    RDouble pressure = primitiveVariables[IP];
    RDouble gasConstant = GetUniversalGasConstant();

    //! To obtain the molecular weight of species.
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();
    //! To obtain the thermodynamic properties of species.
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();
    RDouble massReciprocal = 0.0;// To record the average molecular weight.
    for (int i = 0; i <= nLastSpeciesIndex; ++ i)
    {
        massReciprocal += primitiveVariables[nm + i] / speciesMass[i];
    }

    //! Compute the translation-rotation temperature.
    RDouble transRotationTemperature = (pressure - electronPressure) / (density * gasConstant * massReciprocal);
    return transRotationTemperature;
}

RDouble MixedGas::ComputeTranslationAndRotationTemperature(RDouble *massFraction, RDouble transRotationEnergy)
{
    RDouble transRotationCv = GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(massFraction);
    RDouble formationEnthalpy = GetTotalFormationEnthalpy(massFraction);

    return (transRotationEnergy - formationEnthalpy) / transRotationCv;
}

RDouble MixedGas::ComputeVibrationTemperature(RDouble *massFraction, RDouble vibrationEnergy, RDouble initTemperature/* = 0.0*/)
{
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble e1 = 0.0, e2 = 0.0, T1 = 100.0 / referenceTemperature, T2 = 30000.0 / referenceTemperature;
    RDouble ft1, ft2, Tk = 0.0, Tk1 = 0.0, precision = 1.0e-12;

    int nCount = 0, nMaxStep = 500;
    //! To determine the region where the root exists.
    if (initTemperature <= 0.0)
    {
        e1 = GetMixedGasVibrationEnergy(T1, massFraction);
        e2 = GetMixedGasVibrationEnergy(T2, massFraction);
        if (vibrationEnergy <= e1)
        {
            return T1;
        }
        else if (vibrationEnergy >= e2)
        {
            return T2;
        }

        //! To determine the region where the root exists.
        ft1 = e1 - vibrationEnergy;
        RDouble deltaT = 3000.0 / referenceTemperature;
        while(nCount <= 10)
        {
            nCount ++;
            T2 = nCount * deltaT;
            ft2 = GetMixedGasVibrationEnergy(T2, massFraction) - vibrationEnergy;
            if (fabs(ft2) <= precision)
            {
                return T2;
            }
            else if (ft1 * ft2 < 0.0)
            {
                break;      //! Find the region where the root exists.
            }
            T1 = T2;
            ft1 = ft2;
        }
        if (nCount > 10)     //! Failed to find the region, set the minimum temperature.
        {
            return T2;
        }
        Tk = (T1 + T2) / 2.0;
    }
    else
    {
        nMaxStep = 5;
        Tk = initTemperature;
    }

    //! To compute the vibration temperature using the Newton iteration method.
    nCount = 0;
    RDouble vibrationCv = 0.0, deltaEnergy = 0.0;
    while (nCount < nMaxStep)
    {
        nCount++;
        vibrationCv = GetMixedGasVibrationSpecificHeatAtConstantVolume(Tk, massFraction);
        deltaEnergy = GetMixedGasVibrationEnergy(Tk, massFraction) - vibrationEnergy;
        Tk1 = Tk - deltaEnergy / (vibrationCv + SMALL);
        if (fabs(Tk1 - Tk) <= precision)
        {
            break;
        }
        Tk = Tk1;
    }
    //printf("Number of iteration is %d\n", nCount);
    return Tk1;
}

RDouble MixedGas::ComputeElectronTemperature(RDouble *massFraction, RDouble electronEnergy, RDouble initTemperature)
{
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble e1 = 0.0, e2 = 0.0, T1 = 100.0 / referenceTemperature, T2 = 30000.0 / referenceTemperature;
    RDouble ft1, ft2, Tk = 0.0, Tk1 = 0.0, precision = 1.0e-12;

    //! To determine the region where the root exists.
    int nCount = 0, nMaxStep = 500;
    if (initTemperature <= 0.0)
    {
        e1 = GetMixedGasElectronEnergy(T1, massFraction);
        e2 = GetMixedGasElectronEnergy(T2, massFraction);
        if (electronEnergy <= e1)
        {
            return T1;
        }
        else if (electronEnergy >= e2)
        {
            return T2;
        }

        //! To determine the region where the root exists.
        ft1 = e1 - electronEnergy;
        RDouble deltaT = 3000.0 / referenceTemperature;
        while(nCount <= 10)
        {
            nCount ++;
            T2 = nCount * deltaT;
            ft2 = GetMixedGasElectronEnergy(T2, massFraction) - electronEnergy;
            if (fabs(ft2) <= precision)
            {
                return T2;
            }
            else if (ft1 * ft2 < 0.0)
            {
                break;      //! Find the region where the root exists.
            }
            T1 = T2;
            ft1 = ft2;
        }
        if (nCount > 10)     //! Failed to find the region, set the minimum temperature.
        {
            return T2;
        }
        Tk = (T1 + T2) / 2.0;
    }
    else
    {
        nMaxStep = 5;
        Tk = initTemperature;
    }

    //! To compute the electron temperature using the Newton iteration method.
    nCount = 0;
    RDouble electronCv = 0.0, deltaEnergy = 0.0;
    while (nCount < nMaxStep)
    {
        nCount++;
        electronCv = GetMixedGasElectronSpecificHeatAtConstantVolume(Tk, massFraction);
        deltaEnergy = GetMixedGasElectronEnergy(Tk, massFraction) - electronEnergy;
        Tk1 = Tk - deltaEnergy / (electronCv + SMALL);
        if (fabs(Tk1 - Tk) <= precision)
        {
            break;
        }
        Tk = Tk1;
    }
    //printf("Number of iteration is %d\n", nCount);
    return Tk1;
}

RDouble MixedGas::ComputeVibrationAndElectronTemperature(RDouble *massFraction, RDouble vibrationElectronEnergy, RDouble initTemperature)
{
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble e1 = 0.0, e2 = 0.0, T1 = 100.0 / referenceTemperature, T2 = 30000.0 / referenceTemperature;
    RDouble ft1, ft2, Tk = 0.0, Tk1 = 0.0, precision = 1.0e-12;
    RDouble tvmin = 0.0 / referenceTemperature;

    //! To determine the region where the root exists.
    int nCount = 0, nMaxStep = 500;
    if (initTemperature <= 0.0)
    {
        e1 = GetMixedGasVibrationEnergy(T1, massFraction) + GetMixedGasElectronEnergy(T1, massFraction);
        e2 = GetMixedGasVibrationEnergy(T2, massFraction) + GetMixedGasElectronEnergy(T2, massFraction);

        if (vibrationElectronEnergy <= e1)
        {
            return T1;
        }
        else if (vibrationElectronEnergy >= e2)
        {
            return T2;
        }

        //! To determine the region where the root exists.
        ft1 = e1 - vibrationElectronEnergy;
        RDouble deltaT = 3000.0 / referenceTemperature;
        while(nCount <= 10)
        {
            nCount ++;
            T2 = nCount * deltaT;
            ft2 = GetMixedGasVibrationEnergy(T2, massFraction) + GetMixedGasElectronEnergy(T2, massFraction) - vibrationElectronEnergy;

            if (fabs(ft2) <= precision)
            {
                return T2;
            }
            else if (ft1 * ft2 < 0.0)
            {
                break;     //! Find the region where the root exists.
            }
            T1 = T2;
            ft1 = ft2;
        }
        if (nCount > 10)     //! Failed to find the region, set the minimum temperature.
        {
            return T2;
        }
        Tk = (T1 + T2) / 2.0;
    }
    else
    {
      //  nMaxStep = 100;
        nMaxStep = 5;
        Tk = initTemperature;
    }

    //! To compute the vibration-electron temperature using the Newton iteration method.
    nCount = 0;
    RDouble vibrationElectronCv = 0.0, deltaEnergy = 0.0;

    while (nCount < nMaxStep)
    {
        nCount++;
        vibrationElectronCv = GetMixedGasVibrationSpecificHeatAtConstantVolume(MAX(Tk,tvmin), massFraction);
        vibrationElectronCv += GetMixedGasElectronSpecificHeatAtConstantVolume(MAX(Tk,tvmin), massFraction);

        /*vibrationElectronCv = GetMixedGasVibrationSpecificHeatAtConstantVolume(Tk, massFraction);
        vibrationElectronCv += GetMixedGasElectronSpecificHeatAtConstantVolume(Tk, massFraction);*/

        deltaEnergy = GetMixedGasVibrationEnergy(Tk, massFraction) + GetMixedGasElectronEnergy(Tk, massFraction) - vibrationElectronEnergy;

        RDouble TkLimit = deltaEnergy / (vibrationElectronCv + SMALL);
        if (TkLimit > 0.9 * Tk)
        {
            TkLimit = 0.9 * Tk;
        }
        else if(TkLimit < -10.0 * Tk)
        {
            TkLimit = -10.0 * Tk;
        }
        Tk1 = Tk - TkLimit;

        if (  fabs(TkLimit) <= precision )
        {
            break;
        }
        /*else if ( nCount > 5  && fabs(TkLimit) <= 1.0e-3 )
        {
            break;
        }*/
        
        Tk = Tk1;               
    }

    return Tk1;
}

RDouble MixedGas::ComputeOneTemperatureModelTemperature(RDouble *massFraction, RDouble internalEnergy)
{
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble e1 = 0.0, e2 = 0.0, T1 = 100.0 / referenceTemperature, T2 = 30000.0 / referenceTemperature;
    RDouble ft1, ft2, precision = 1.0e-12;
    e1 = GetMixedGasVibrationEnergy(T1, massFraction) + GetMixedGasElectronEnergy(T1, massFraction);
    e2 = GetMixedGasVibrationEnergy(T2, massFraction) + GetMixedGasElectronEnergy(T2, massFraction);
    e1 += GetMixedGasTranslationAndRotationEnergy(T1, massFraction);
    e2 += GetMixedGasTranslationAndRotationEnergy(T2, massFraction);
    if (internalEnergy <= e1)
    {
        return T1;
    }
    else if (internalEnergy >= e2)
    {
        return T2;
    }

    //! To determine the region where the root exists.
    ft1 = e1 - internalEnergy;
    RDouble deltaT = 3000.0 / referenceTemperature;
    int nCount = 0, nMaxStep = 500;
    while(nCount <= 10)
    {
        nCount ++;
        T2 = nCount * deltaT;
        ft2 = GetMixedGasVibrationEnergy(T2, massFraction) + GetMixedGasElectronEnergy(T2, massFraction);
        ft2 += GetMixedGasTranslationAndRotationEnergy(T2, massFraction) - internalEnergy;
        if (fabs(ft2) <= precision)
        {
            return T2;
        }
        else if (ft1 * ft2 < 0.0)
        {
            break;      //! Find the region where the root exists.
        }
        T1 = T2;
        ft1 = ft2;
    }
    if (nCount > 10)     //! Failed to find the region, set the minimum temperature.
    {
        return T2;
    }

    //! To compute the vibration-electron temperature using the Newton iteration method.
    nCount = 0;
    RDouble Tk = (T1 + T2) / 2.0, Tk1 = 0.0, gasMixtureCv = 0.0, deltaEnergy = 0.0;
    while (nCount < nMaxStep)
    {
        nCount++;
        gasMixtureCv = GetMixedGasVibrationSpecificHeatAtConstantVolume(Tk, massFraction);
        gasMixtureCv += GetMixedGasElectronSpecificHeatAtConstantVolume(Tk, massFraction);
        gasMixtureCv += GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(massFraction);

        deltaEnergy = GetMixedGasVibrationEnergy(Tk, massFraction) + GetMixedGasElectronEnergy(Tk, massFraction);
        deltaEnergy += GetMixedGasTranslationAndRotationEnergy(Tk, massFraction) - internalEnergy;
        Tk1 = Tk - deltaEnergy / gasMixtureCv;
        if (fabs(Tk1 - Tk) <= precision)
        {
            break;
        }
        Tk = Tk1;
    }
    //printf("Number of iteration is %d\n", nCount);
    return Tk1;
}

RDouble MixedGas::ComputeNonequilibrumTemperatureViaBisectionMethod(RDouble *massFraction, RDouble internalEnergy)
{
    RDouble gasConstant = this->GetUniversalGasConstant();
    RDouble mass = this->GetMixedGasMolecularWeight(massFraction);

    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble gasR = gasConstant / mass;
    RDouble Tm, T1, T2, e1, e2, ftm, ft1, ft2, precision = 1.0e-12;
    T1 = 100.0 / referenceTemperature;
    T2 = 30000.0 / referenceTemperature;
    e1 = GetMixtureGasEnthalpy(massFraction, T1) - gasR * T1;
    e2 = GetMixtureGasEnthalpy(massFraction, T2) - gasR * T2;
    if (internalEnergy <= e1)
    {
        return T1;
    }
    else if (internalEnergy >= e2)
    {
        return T2;
    }

    //! To determine the region where the root exists.
    ft1 = e1 - internalEnergy;
    RDouble deltaT = 3000.0 / referenceTemperature;
    int nCount = 0, nMaxStep = 500;
    while(nCount <= 10)
    {
        nCount ++;
        T2 = nCount * deltaT;
        ft2 = GetMixtureGasEnthalpy(massFraction, T2) - gasR * T2 - internalEnergy;
        if (fabs(ft2) <= precision)
        {
            return T2;
        }
        else if (ft1 * ft2 < 0.0)
        {
            break;      //! Find the region where the root exists.
        }
        T1 = T2;
        ft1 = ft2;
    }
    if (nCount > 10)     //! Failed to find the region, set the minimum temperature.
    {
        return T2;
    }

    //! To compute the vibration temperature using the bisection method.
    nCount = 0;
    while (nCount < nMaxStep)
    {
        nCount++;
        if (nCount >= nMaxStep)
        {
            break;
        }

        Tm = 0.5 * (T1 + T2);
        ftm = GetMixtureGasEnthalpy(massFraction, Tm) - gasR * Tm - internalEnergy;
        if (fabs(ftm) <= precision)     //! The root is found.
        {
            break;
        }

        if (ft1 * ftm < 0.0)
        {
            T2 = Tm;
            ft2 = ftm;
        }
        if (ft2 * ftm < 0.0)
        {
            T1 = Tm;
            ft1 = ftm;
        }
        if (fabs(T1 - T2) <= precision)    //! The root approaches to T1 or T2.
        {
            Tm = T1;
            break;
        }
    }

    if (nCount >= nMaxStep)
    {
        Tm = 0.5 * (T1 + T2);
    }
    return Tm;
}

RDouble MixedGas::ComputeNonequilibrumTemperature(RDouble *massFraction, RDouble internalEnergy, RDouble initTemperature)
{
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble e1 = 0.0, e2 = 0.0, T1 = 100.0 / referenceTemperature, T2 = 30000.0 / referenceTemperature;
    RDouble ft1, ft2, Tk = 0.0, Tk1 = 0.0, precision = 1.0e-12;
    RDouble gasConstant = this->GetUniversalGasConstant();
    RDouble mass = this->GetMixedGasMolecularWeight(massFraction);
    RDouble gasR = gasConstant / mass;

    //! To determine the region where the root exists.
    int nCount = 0, nMaxStep = 500;
    if (initTemperature <= 0.0)
    {
        e1 = GetMixtureGasEnthalpy(massFraction, T1) - gasR * T1;
        e2 = GetMixtureGasEnthalpy(massFraction, T2) - gasR * T2;
        if (internalEnergy <= e1)
        {
            return T1;
        }
        else if (internalEnergy >= e2)
        {
            return T2;
        }

        //! To determine the region where the root exists.
        ft1 = e1 - internalEnergy;
        RDouble deltaT = 3000.0 / referenceTemperature;
        while(nCount <= 10)
        {
            nCount ++;
            T2 = nCount * deltaT;
            ft2 = GetMixtureGasEnthalpy(massFraction, T2) - gasR * T2 - internalEnergy;
            if (fabs(ft2) <= precision)
            {
                return T2;
            }
            else if (ft1 * ft2 <= 0.0)
            {
                break;     //! Find the region where the root exists.
            }
            T1 = T2;
            ft1 = ft2;
        }
        if (nCount > 10)     //! Failed to find the region, set the minimum temperature.
        {
            return T2;
        }
        Tk = (T1 + T2) / 2.0;
    }
    else
    {
        nMaxStep = 5;
        Tk = initTemperature;
    }

    //! To compute the vibration-electron temperature using the Newton iteration method.
    nCount = 0;
    RDouble gasMixtureCv = 0.0, deltaEnergy = 0.0;
    while (nCount < nMaxStep)
    {
        nCount++;
        ComputeConstantVolumeSpecificHeatByMassFraction(massFraction, Tk, gasMixtureCv);
        deltaEnergy = GetMixtureGasEnthalpy(massFraction, Tk) - gasR * Tk - internalEnergy;
        Tk1 = Tk - deltaEnergy / gasMixtureCv;
        if (fabs(Tk1 - Tk) <= precision)
        {
            break;
        }
        Tk = Tk1;
    }
    //printf("Number of iteration is %d\n", nCount);
    return Tk1;
}

RDouble MixedGas::GetElectronPressure(RDouble *primitiveVariables, RDouble electronTemperature)
{
    RDouble density = primitiveVariables[IDX::IR];
    RDouble gasConstant = this->GetUniversalGasConstant();
    //RDouble referenceMass = referenceParameterDimensional.GetAverageMolecularWeight();
    //! To obtain the molecular weight of species.
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();

    //! Search the index of the electron.
    RDouble electronMassFraction = 0.0;
    //RDouble electronMass = 5.486e-7 / referenceMass;
    if (nElectronIndex >= 0)
    {
        electronMassFraction = primitiveVariables[nm + nElectronIndex];
        //electronMass = speciesMass[nElectronIndex];
    }

    //! Compute electron pressure.
    RDouble electronPressure = (density * gasConstant * electronTemperature) * electronMassFraction / electronMass;
    return electronPressure;
}

void MixedGas::ComputeSpeciesTranslationAndRotationHeatConductivity(RDouble *viscosity, RDouble *conductivity)
{
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    RDouble gasConstant = GetUniversalGasConstant();
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();

    RDouble *cvts = new RDouble[this->numberOfSpecies];
    RDouble *cvrs = new RDouble[this->numberOfSpecies];

    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        //! Compute translation specific heat at constant volume.
        cvts[i] = 1.5 * gasConstant / speciesMass[i];

        //! Compute rotation specific heat at constant volume.
        if (thermodynamicProperties[i].nType >= 2)     //! Polyatom species.
        {
            cvrs[i] = gasConstant / speciesMass[i]; 
        }
        else     //! Monoatom species or electron.
        {
            cvrs[i] = 0.0;
        }

        //! Compute translation-rotation heat conductivity of species.
        conductivity[i] = viscosity[i] * (2.5 * cvts[i] + cvrs[i]);
    }

    delete []cvts;
    delete []cvrs;
}

void MixedGas::ComputeSpeciesVibrationHeatConductivity(RDouble *viscosity, RDouble *speciesVibrationCv, RDouble *conductivity)
{
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        conductivity[i] = viscosity[i] * speciesVibrationCv[i];
    }
}

void MixedGas::ComputeSpeciesElectronHeatConductivity(RDouble *viscosity, RDouble *speciesElectronCv, RDouble *conductivity)
{
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        conductivity[i] = viscosity[i] * speciesElectronCv[i];
    }
}

RDouble MixedGas::GetMixedGasTranslationAndRotationHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature)
{
    RDouble *viscosity = new RDouble[this->numberOfSpecies];
    RDouble *conductivity = new RDouble[this->numberOfSpecies];
    RDouble *massFraction = new RDouble[this->numberOfSpecies];
    RDouble *moleFraction = new RDouble[this->numberOfSpecies];
    RDouble *phi = new RDouble[this->numberOfSpecies];

    //! Obtain the molecular weight of each species.
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();
    for (int s = 0; s < this->numberOfSpecies; ++ s)
    {
        massFraction[s] = primitiveVariables[nm + s];
    }
    //! Obtain the mole fraction of each species.
    this->ComputeMoleFractionByMassFraction(massFraction, moleFraction);

    //! Compute viscosity of species.
    int nType = 0;      //! 0 for Blottnner curve fit method, 1 for Lennard-Johns method.
    if (nType == 0)
    {
        ComputeSpeciesViscosityWithCurveFitMethod(transRotationTemperature, electronTemperature, viscosity);
    }
    else
    {
        RDouble electronPressure = primitiveVariables[IDX::IP];
        RDouble gasConstant = GetUniversalGasConstant();
        if (nElectronIndex >= 0)
        {
            electronPressure = primitiveVariables[IDX::IR] * electronTemperature * gasConstant * massFraction[nElectronIndex] / speciesMass[nElectronIndex];
        }

        ComputeSpeciesViscosityWithLennardJonesMethod(transRotationTemperature, electronPressure, viscosity);
    }

    //! Compute the translation-rotation heat conductivity.
    ComputeSpeciesTranslationAndRotationHeatConductivity(viscosity, conductivity);

    //! Obtain the partition function.
    GetPartitionFunctionValues(viscosity, speciesMass, moleFraction, phi);

    RDouble lambda = GetMixedGasHeatConductivityWithWilkeFormula(moleFraction, conductivity, phi);

    delete []viscosity;
    delete []conductivity;
    delete []massFraction;
    delete []moleFraction;
    delete []phi;
    //! Return the heat conductivity of translation and rotation.
    return lambda;
}

RDouble MixedGas::GetMixedGasVibrationHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature)
{
    RDouble *viscosity = new RDouble[this->numberOfSpecies];
    RDouble *speciesVibrationCv = new RDouble[this->numberOfSpecies];
    RDouble *conductivity = new RDouble[this->numberOfSpecies];
    RDouble *massFraction = new RDouble[this->numberOfSpecies];
    RDouble *moleFraction = new RDouble[this->numberOfSpecies];
    RDouble *phi = new RDouble[this->numberOfSpecies];

    //! Obtain the molecular weight of each species.
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        massFraction[i] = primitiveVariables[nm + i];
    }
    //! Obtain the mole fraction of each species.
    this->ComputeMoleFractionByMassFraction(massFraction, moleFraction);

    //! Compute viscosity of species.
    int nType = 0;     //! 0 for Blottnner curve fit method, 1 for Lennard-Johns method.
    if (nType == 0)
    {
        ComputeSpeciesViscosityWithCurveFitMethod(transRotationTemperature, electronTemperature, viscosity);
    }
    else
    {
        RDouble electronPressure = primitiveVariables[IDX::IP];
        RDouble gasConstant = GetUniversalGasConstant();
        if (nElectronIndex >= 0)
        {
            electronPressure = primitiveVariables[IDX::IR] * electronTemperature * gasConstant * massFraction[nElectronIndex] / speciesMass[nElectronIndex];
        }

        ComputeSpeciesViscosityWithLennardJonesMethod(transRotationTemperature, electronPressure, viscosity);
    }

    //! Compute vibration specific heat at constant volume.
    ComputeVibrationSpecificHeatAtConstantVolume(vibrationTemperature, speciesVibrationCv);

    //! Compute the vibration heat conductivity.
    ComputeSpeciesVibrationHeatConductivity(viscosity, speciesVibrationCv, conductivity);

    //! Obtain the partition function.
    GetPartitionFunctionValues(viscosity, speciesMass, moleFraction, phi);

    RDouble lambda = GetMixedGasHeatConductivityWithWilkeFormula(moleFraction, conductivity, phi);

    delete []viscosity;
    delete []speciesVibrationCv;
    delete []conductivity;
    delete []massFraction;
    delete []moleFraction;
    delete []phi;
    //! Return the heat conductivity of vibration.
    return lambda;
}

RDouble MixedGas::GetMixedGasElectronHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature)
{
    RDouble *viscosity = new RDouble[this->numberOfSpecies];
    RDouble *speciesElectronCv = new RDouble[this->numberOfSpecies];
    RDouble *conductivity = new RDouble[this->numberOfSpecies];
    RDouble *massFraction = new RDouble[this->numberOfSpecies];
    RDouble *moleFraction = new RDouble[this->numberOfSpecies];
    RDouble *phi = new RDouble[this->numberOfSpecies];

    //! Obtain the molecular weight of each species.
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        massFraction[i] = primitiveVariables[nm + i];
    }
    //! Obtain the mole fraction of each species.
    this->ComputeMoleFractionByMassFraction(massFraction, moleFraction);

    //! Compute viscosity of species.
    int nType = 0;      //! 0 for Blottnner curve fit method, 1 for Lennard-Johns method.
    if (nType == 0)
    {
        ComputeSpeciesViscosityWithCurveFitMethod(transRotationTemperature, electronTemperature, viscosity);
    }
    else
    {
        RDouble electronPressure = primitiveVariables[IDX::IP];
        RDouble gasConstant = GetUniversalGasConstant();
        if (nElectronIndex >= 0)
        {
            electronPressure = primitiveVariables[IDX::IR] * electronTemperature * gasConstant * massFraction[nElectronIndex] / speciesMass[nElectronIndex];
        }

        ComputeSpeciesViscosityWithLennardJonesMethod(transRotationTemperature, electronPressure, viscosity);
    }

    //! Compute vibration specific heat at constant volume.
    ComputeElectronSpecificHeatAtConstantVolume(electronTemperature, speciesElectronCv);

    //! Compute the electron heat conductivity.
    ComputeSpeciesElectronHeatConductivity(viscosity, speciesElectronCv, conductivity);

    //! Obtain the partition function.
    GetPartitionFunctionValues(viscosity, speciesMass, moleFraction, phi);

    RDouble lambda = GetMixedGasHeatConductivityWithWilkeFormula(moleFraction, conductivity, phi);

    delete []viscosity;
    delete []speciesElectronCv;
    delete []conductivity;
    delete []massFraction;
    delete []moleFraction;
    delete []phi;
    //! Return the heat conductivity of vibration.
    return lambda;
}

void MixedGas::ComputeMixedGasTranslationVibrationElectronHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble &transRotationConductivity, RDouble &vibrationConductivity, RDouble &electronConductivity)
{
    RDouble *viscosity = new RDouble[this->numberOfSpecies];
    RDouble *conductivity = new RDouble[this->numberOfSpecies];
    RDouble *massFraction = new RDouble[this->numberOfSpecies];
    RDouble *moleFraction = new RDouble[this->numberOfSpecies];
    RDouble *phi = new RDouble[this->numberOfSpecies];
    RDouble *speciesVibrationCv = new RDouble[this->numberOfSpecies];
    RDouble *speciesElectronCv = new RDouble[this->numberOfSpecies];

    //! Obtain the molecular weight of each species.
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();
    for (int s = 0; s < this->numberOfSpecies; ++ s)
    {
        massFraction[s] = primitiveVariables[nm + s];
    }
    //! Obtain the mole fraction of each species.
    this->ComputeMoleFractionByMassFraction(massFraction, moleFraction);

    //! Compute viscosity of species.
    int nType = 0;     //! 0 for Blottnner curve fit method, 1 for Lennard-Johns method.
    if (nType == 0)
    {
        ComputeSpeciesViscosityWithCurveFitMethod(transRotationTemperature, electronTemperature, viscosity);
    }
    else
    {
        RDouble electronPressure = primitiveVariables[IDX::IP];
        RDouble gasConstant = GetUniversalGasConstant();
        if (nElectronIndex >= 0)
        {
            electronPressure = primitiveVariables[IDX::IR] * electronTemperature * gasConstant * massFraction[nElectronIndex] / speciesMass[nElectronIndex];
        }

        ComputeSpeciesViscosityWithLennardJonesMethod(transRotationTemperature, electronPressure, viscosity);
    }

    //! Compute the translation-rotation heat conductivity.
    ComputeSpeciesTranslationAndRotationHeatConductivity(viscosity, conductivity);

    //! Obtain the partition function.
    GetPartitionFunctionValues(viscosity, speciesMass, moleFraction, phi);

    transRotationConductivity = GetMixedGasHeatConductivityWithWilkeFormula(moleFraction, conductivity, phi);

    //! Compute vibration specific heat at constant volume.
    ComputeVibrationSpecificHeatAtConstantVolume(vibrationTemperature, speciesVibrationCv);

    //! Compute the vibration heat conductivity.
    ComputeSpeciesVibrationHeatConductivity(viscosity, speciesVibrationCv, conductivity);

    vibrationConductivity = GetMixedGasHeatConductivityWithWilkeFormula(moleFraction, conductivity, phi);

    //! Compute vibration specific heat at constant volume.
    ComputeElectronSpecificHeatAtConstantVolume(electronTemperature, speciesElectronCv);

    //! Compute the electron heat conductivity.
    ComputeSpeciesElectronHeatConductivity(viscosity, speciesElectronCv, conductivity);

    electronConductivity = GetMixedGasHeatConductivityWithWilkeFormula(moleFraction, conductivity, phi);
    
    delete []viscosity;
    delete []conductivity;
    delete []massFraction;
    delete []moleFraction;
    delete []phi;
    delete []speciesVibrationCv;
    delete []speciesElectronCv;
}

RDouble MixedGas::GetMixedGasHeatConductivityWithWilkeFormula(RDouble *moleFraction, RDouble *conductivity, RDouble *phi)
{
    RDouble lambda = 0.0;
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        lambda += moleFraction[i] * conductivity[i] / phi[i];
    }
    return lambda;
}

void MixedGas::ComputeSpeciesViscosityWithCurveFitMethod(RDouble transRotationTemperature, RDouble electronTemperature, RDouble *speciesViscosity)
{
    //! Obtain the reference viscosity.
    RDouble referenceViscosity = this->referenceParameterDimensional.GetViscosity();

    //! Obtain the dimensional value of viscosity.
    ComputeSpeciesDimensionalViscosityWithCurveFitMethod(transRotationTemperature, electronTemperature, speciesViscosity);

    //! Translate to the non-dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesViscosity[i] /= referenceViscosity;
    }
}

void MixedGas::ComputeSpeciesDimensionalViscosityWithCurveFitMethod(RDouble transRotationTemperature, RDouble electronTemperature, RDouble *speciesViscosity)
{
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble dimensionalTransRotationTemperature, dimensionalElectronTemperature;
    dimensionalTransRotationTemperature = transRotationTemperature * referenceTemperature;
    dimensionalElectronTemperature = electronTemperature * referenceTemperature;
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();

    //! Obtain the coefficients of the curve fitting.
    RDouble *coefA = CurveFits->GetBlotterCurveFittingCoefficientA();
    RDouble *coefB = CurveFits->GetBlotterCurveFittingCoefficientB();
    RDouble *coefC = CurveFits->GetBlotterCurveFittingCoefficientC();
    RDouble *coefD = CurveFits->GetBlotterCurveFittingCoefficientD();
    RDouble *coefE = CurveFits->GetBlotterCurveFittingCoefficientE();

    RDouble lnt1, lnt2, lnt3, lnt4, tmp;
    lnt1 = log(dimensionalTransRotationTemperature);     //! Polyatom and monoatom.
    lnt2 = lnt1 * lnt1;
    lnt3 = lnt1 * lnt2;
    lnt4 = lnt1 * lnt3;
    for (int i = 0; i <= nLastSpeciesIndex; ++ i)
    {
        //! To compute viscosity using the curve fitting method.
        tmp = coefA[i] * lnt4 + coefB[i] * lnt3 + coefC[i] * lnt2 + coefD[i] * lnt1 + coefE[i];

        //! Compute the dimensional value of viscosity.
        speciesViscosity[i] = 0.1 * exp(tmp);
    }
    if (nElectronIndex >= 0)// electron.
    {
        lnt1 = log(dimensionalElectronTemperature);
        lnt2 = lnt1 * lnt1;
        lnt3 = lnt1 * lnt2;
        lnt4 = lnt1 * lnt3;
        tmp = coefA[nElectronIndex] * lnt4 + coefB[nElectronIndex] * lnt3 + coefC[nElectronIndex] * lnt2 + coefD[nElectronIndex] * lnt1 + coefE[nElectronIndex];
        speciesViscosity[nElectronIndex] = 0.1 * exp(tmp);
    }
}

void MixedGas::ComputeSpeciesViscosityWithLennardJonesMethod(RDouble transRotationTemperature, RDouble electronPressure, RDouble *speciesViscosity)
{
    //! Obtain the reference viscosity.
    RDouble referenceViscosity = this->referenceParameterDimensional.GetViscosity();

    //! Obtain the dimensional value of viscosity.
    ComputeSpeciesDimensionalViscosityWithLennardJonesMethod(transRotationTemperature, electronPressure, speciesViscosity);

    //! Translate to the non-dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesViscosity[i] /= referenceViscosity;
    }
}

void MixedGas::ComputeSpeciesDimensionalViscosityWithLennardJonesMethod(RDouble transRotationTemperature, RDouble electronPressure, RDouble *speciesViscosity)
{
    //! Obtain the reference temperature.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble dimensionalTransRotationTemperature = transRotationTemperature * referenceTemperature;

    //! The average area of collision between species.
    RDouble *collisionArea = new RDouble[this->numberOfSpecies];
    ComputeAverageCollisionAreaOmega22(transRotationTemperature, electronPressure, collisionArea);

    //! Obtain the molecular weights of species.
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();
    //! Translate to the dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesViscosity[i] = 8.4411e-8 * PI * sqrt(speciesMass[i] * dimensionalTransRotationTemperature) / collisionArea[i];
    }
    delete []collisionArea;
}

void MixedGas::ComputeAverageCollisionAreaOmega22(RDouble transRotationTemperature, RDouble electronPressure, RDouble *piOmega)
{
    //! Obtain the reference temperature and reference pressure.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referencePressure = this->referenceParameterDimensional.GetPressure();
    RDouble dimensionalTransRotationTemperature, temperature, lnt, tmp;
    dimensionalTransRotationTemperature = transRotationTemperature * referenceTemperature;
    temperature = dimensionalTransRotationTemperature / 1000.0;
    lnt = log(dimensionalTransRotationTemperature);

    RDouble p0, pe, pem, pressure, lnpe;
    pe  = electronPressure * referencePressure;
    p0  = 101325.0; //Pa
    pem = p0 * 0.0975 * pow(temperature, 4);
    pressure = p0 / pe;

    int *ionType = this->molecularProperty->GetIonTypeOfSpecies();
    RDouble *coefAss, *coefBss, *coefCss, *coefDss;
    coefAss = this->CurveFits->GetOmegaCurveFittingCoefficientA();
    coefBss = this->CurveFits->GetOmegaCurveFittingCoefficientB();
    coefCss = this->CurveFits->GetOmegaCurveFittingCoefficientC();
    coefDss = this->CurveFits->GetOmegaCurveFittingCoefficientD();

    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        if (ionType[i] == 0)     //! The atom or molecular.
        {
            lnpe = 1.0;
        }
        else
        {
            if (fabs(pe - pem) <= 1.0e-8)
            {
                lnpe = 1.0;
            }
            else     //! Modify the coefficient.
            {
                lnpe = 0.5 * log(2.09e-2 * pow(temperature, 4) * pressure + 1.52 * pow(temperature, 8.0 / 3.0) * pow(pressure, 2.0 / 3.0));
            }
        }
        tmp = coefAss[i] * pow(lnt, 3) + coefBss[i] * pow(lnt, 2) + coefCss[i] * lnt + coefDss[i];
        piOmega[i] = 1.0e-20 * lnpe * exp(tmp);
    }
}

RDouble MixedGas::GetMixedGasViscosityWithWilkeFormula(RDouble *moleFraction, RDouble *viscosity, RDouble *phi)
{
    RDouble mixedViscosity = 0.0;
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        mixedViscosity += moleFraction[i] * viscosity[i] / phi[i];
    }
    return mixedViscosity;
}

void MixedGas::ComputeSpeciesMassDiffusionCoefficient(RDouble *massFraction, RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble *speciesDiffusionCoef)
{
    //! The function returns the diffusion coefficients of each species which saved in rDs = (1 - Cs)/(1- Xs) * (mul/Scl + mut/Sct).
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    int *ionType = this->molecularProperty->GetIonTypeOfSpecies();

    //! To compute the mole fractions of species.
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();
    this->ComputeMoleFractionByMassFraction(massFraction, moleFractions);

    //! Obtain the Schmidt number.
    RDouble *oSchmidtLaminar = GetLaminarSchmidtNumberReciprocal();
    RDouble *oSchmidtTurbulence = GetTurbulentSchmidtNumberReciprocal();

    RDouble tmp1 = 0.0, tmp2 = 0.0;
    int nIndex = -1;
    for(int i = 0; i < this->numberOfSpecies; ++ i)
    {
        if (thermodynamicProperties[i].nType == 0)
        {
            nIndex = i;      //! To save the index of the electron.
            continue;      //! Ignore the electron.
        }

        if (1.0 - moleFractions[i] <= 1.0e-20)
        {
            speciesDiffusionCoef[i] = 0.0;
        }
        else
        {
            speciesDiffusionCoef[i] = (1.0 - massFraction[i]) / (1.0 - moleFractions[i]) * (viscosityLaminar * oSchmidtLaminar[i] + viscosityTurbulence * oSchmidtTurbulence[i]);
        }

        if (ionType[i] == 0)
        {
            continue;     //! Ignore the atom and molecular species.
        }
        else
        {
            tmp1 += speciesDiffusionCoef[i] * massFraction[i] / speciesMass[i];
            tmp2 += massFraction[i] / speciesMass[i];
        }
    }

    //! Compute the mass diffusion coefficient of electron.
    if (nIndex >= 0)
    {
        speciesDiffusionCoef[nIndex] = tmp1 / MAX(tmp2, 1.0e-20);
    }
}

RDouble MixedGas::GetMassDiffusionCoefficientWithCLNModel(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature)
{
    RDouble *massFraction = new RDouble[this->numberOfSpecies]();
    //! Obtain the mass fractions of species.
    for(int i = 0; i < this->numberOfSpecies; ++ i)
    {
        massFraction[i] = primitiveVariables[nm + i];
    }

    RDouble gasConstant = this->GetUniversalGasConstant();
    RDouble mixedMass = this->GetMixedGasMolecularWeight(massFraction);
    RDouble transRotationCv = GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(massFraction);
    RDouble transRotationCp = transRotationCv + gasConstant / mixedMass;
    RDouble lewisNumber = 1.4;     //! The default value of Lewis number.
    RDouble mixedConductivity = GetMixedGasTranslationAndRotationHeatConductivity(primitiveVariables, transRotationTemperature, electronTemperature);
    delete []massFraction;    massFraction = NULL;

    //! Return the equivalent mass diffusion coefficients.
    return (mixedConductivity * lewisNumber / transRotationCp);
}

void MixedGas::GetElectronPressurePartialDerivatives(RDouble *primitiveVars, RDouble electronTemperature, RDouble &gama, RDouble *speciesGama)
{
    if (nElectronIndex < 0)
    {
        gama = 0.0;
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            speciesGama[s] = 0.0;
        }
        return;
    }

    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        //! Initialization.
        massFractions[s] = primitiveVars[nm + s];
    }

    //! Compute the specific heat at constant volume.
    RDouble mixGasCv = 0.0;
    if (ntmodel == 2)
    {
        mixGasCv = GetMixedGasVibrationSpecificHeatAtConstantVolume(electronTemperature, massFractions);
        mixGasCv += GetMixedGasElectronSpecificHeatAtConstantVolume(electronTemperature, massFractions);
    }
    else if (ntmodel == 3)
    {
        mixGasCv = GetMixedGasElectronSpecificHeatAtConstantVolume(electronTemperature, massFractions);
    }
    RDouble gasConstant = GetUniversalGasConstant();
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocal();

    //! Compute partial derivative of electron pressure to the vibration or electron energy.
    gama = gasConstant * massFractions[nElectronIndex] * speciesMassReciprocal[nElectronIndex] / (mixGasCv + SMALL);

    //! Compute partial derivatives of electron pressure to the species density.
    for (int s = 0; s < nLastSpeciesIndex; ++ s)
    {
        speciesGama[s] = gasConstant * electronTemperature * ionType[s] * speciesMassReciprocal[s];
    }

}

void MixedGas::GetPressurePartialDerivatives(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble electronTemperature,
                                             RDouble &alpha, RDouble &beta, RDouble *speciesBeta)
{
    //! Initialization.
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        massFractions[s] = primitiveVars[nm + s];
    }
    RDouble gasConstant = GetUniversalGasConstant();
    RDouble mixGasMass = this->GetMixedGasMolecularWeight(massFractions);
    RDouble transRotationCv = GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(massFractions);
    //! Compute the partial derivative of pressure to internal energy.
    alpha = gasConstant / (mixGasMass * transRotationCv);

    //! Compute translation-rotation energies of species.
    ComputeTranslationAndRotationEnergy(transRotationTemperature, speciesEnergyTr);

    //! Compute the energy terms.
    using namespace IDX;
    RDouble uVelocity = primitiveVars[IU];
    RDouble vVelocity = primitiveVars[IV];
    RDouble wVelocity = primitiveVars[IW];
    RDouble squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;

    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocal();

    //! Compute the coefficient called beta.
    RDouble energy = 0.5 * squareVelocity - speciesEnergyTr[nLastSpeciesIndex];
    beta = alpha * energy + gasConstant * transRotationTemperature * speciesMassReciprocal[nLastSpeciesIndex];

    //! Compute the other coefficients saved in the array of speciesBeta.
    RDouble coefA = 0.0, speciesEnergy = 0.0;
    speciesBeta[nLastSpeciesIndex] = 0.0;
    for (int s = 0; s < nLastSpeciesIndex; ++ s)
    {
        coefA = ionType[s] * speciesMassReciprocal[s] * electronMass;
        speciesEnergy = (1.0 + coefA) * speciesEnergyTr[nLastSpeciesIndex] - speciesEnergyTr[s];

        speciesBeta[s]  = alpha * speciesEnergy;
        speciesBeta[s] += gasConstant * transRotationTemperature * (speciesMassReciprocal[s] - speciesMassReciprocal[nLastSpeciesIndex]);
        speciesBeta[s] += coefA * gasConstant * (electronTemperature / electronMass - transRotationTemperature * speciesMassReciprocal[nLastSpeciesIndex]);
    }

}

void MixedGas::GetMultiTemperatureModelPartialDerivatives(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble electronTemperature,
                                                             RDouble &alpha, RDouble &beta, RDouble &gama, RDouble *speciesBeta, RDouble *speciesGama)
{
    //! Compute partial derivatives of electron pressure to conservative variables.
    this->GetElectronPressurePartialDerivatives(primitiveVars, electronTemperature, gama, speciesGama);
    //! Compute partial derivatives of pressure to conservative variables.
    this->GetPressurePartialDerivatives(primitiveVars, transRotationTemperature, electronTemperature, alpha, beta, speciesBeta);
}

void MixedGas::GetMultiTemperatureModelPartialDerivativesR(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble electronTemperature,
                                                          RDouble &alpha, RDouble &beta, RDouble &gama, RDouble *speciesBeta, RDouble *speciesGama,
                                                          RDouble squareVelocity, RDouble *speciesEtrs, RDouble totalCvtr,RDouble totalCvv, RDouble totalCve)
{
    //! Compute partial derivatives of electron pressure to conservative variables.
    //this->GetElectronPressurePartialDerivatives(primitiveVars, electronTemperature, gama, speciesGama);

    RDouble gasConstant = GetUniversalGasConstant();
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    RDouble *speciesMassReciprocal = molecularProperty->GetMolecularWeightReciprocal();

    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        //! Initialization.
        massFractions[s] = primitiveVars[nm + s];
    }

    if (nElectronIndex < 0)
    {
        gama = 0.0;
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            speciesGama[s] = 0.0;
        }
    }    
    else
    {
        RDouble mixGasCv = 0.0;
        if (ntmodel == 2)
        {
            mixGasCv = totalCvv + totalCve;
        }
        else if (ntmodel == 3)
        {
            mixGasCv = totalCve;
        }
        //! Compute partial derivative of electron pressure to the vibration or electron energy.
        gama = gasConstant * massFractions[nElectronIndex] * speciesMassReciprocal[nElectronIndex] / (mixGasCv + SMALL);

        //! Compute partial derivatives of electron pressure to the species density.
        for (int s = 0; s < nLastSpeciesIndex; ++ s)
        {
            speciesGama[s] = gasConstant * electronTemperature * ionType[s] * speciesMassReciprocal[s];
        }
    }

    RDouble mixGasMass = 0.0;
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        mixGasMass += massFractions[i] * speciesMassReciprocal[i];
    }

    //! Compute the partial derivative of pressure to internal energy.
    alpha = gasConstant * mixGasMass / totalCvtr; 
    beta = alpha * (0.5 * squareVelocity - speciesEtrs[nLastSpeciesIndex]) + gasConstant * transRotationTemperature * speciesMassReciprocal[nLastSpeciesIndex];

    //! Compute the other coefficients saved in the array of speciesBeta.
    RDouble coefA, speciesEnergy;
    for (int s = 0; s < nLastSpeciesIndex; ++ s)
    {
        coefA = ionType[s] * speciesMassReciprocal[s] * electronMass;
        speciesEnergy = (1.0 + coefA) * speciesEtrs[nLastSpeciesIndex] - speciesEtrs[s];

        speciesBeta[s]  = alpha * speciesEnergy;
        speciesBeta[s] += gasConstant * transRotationTemperature * (speciesMassReciprocal[s] - speciesMassReciprocal[nLastSpeciesIndex]);
        speciesBeta[s] += coefA * gasConstant * (electronTemperature / electronMass - transRotationTemperature * speciesMassReciprocal[nLastSpeciesIndex]);
    }
    speciesBeta[nLastSpeciesIndex] = 0.0;

}


RDouble MixedGas::ComputePressureDerivativeToTotalDensity(RDouble *primitiveVars)
{
    RDouble *massFraction = new RDouble[numberOfSpecies];
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        massFraction[s] = primitiveVars[nm + s];
    }
    RDouble gasConstant = GetUniversalGasConstant();
    RDouble mixGasMass = this->GetMixedGasMolecularWeight(massFraction);
    RDouble transRotationCv = GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(massFraction);
    //! Compute the partial derivative of pressure to internal energy.
    RDouble alpha = gasConstant / (mixGasMass * transRotationCv);

    delete []massFraction;
    return alpha;
}

RDouble MixedGas::GetMixedGasEnthalpy(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature)
{
    int nChemical = GetChemicalType();
    RDouble staticEnthaly = 0.0;
    if (nChemical == 0)
    {
        staticEnthaly = 3.5 * primitiveVars[IDX::IP] / primitiveVars[IDX::IR];
    }
    else
    {
        RDouble temperatures[3] = {transRotationTemperature, vibrationTemperature, electronTemperature};
        staticEnthaly = ComputeMixedGasEnthalpy(primitiveVars, temperatures);
    }

    //delete []massFraction;
    return staticEnthaly;
}

RDouble MixedGas::ComputeMolecularWeightReciprocal(RDouble *prim)
{
    RDouble omav = 1.0;
    ComputeMolecularWeightReciprocal(prim, omav);
    return omav;
}

RDouble MixedGas::ComputeMolecularWeight(RDouble *prim)
{
    RDouble omav = 1.0;
    ComputeMolecularWeightReciprocal(prim, omav);
    return (1.0 / omav);
}

void MixedGas::ComputeSpeciesMassDiffusionFlux(RDouble *primitiveVars, RDouble *dfsdx, RDouble *dfsdy, RDouble *dfsdz, RDouble transRotationTemperature, RDouble vibrationTemperatuare, RDouble electronTemperature,
                                     RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble nx, RDouble ny, RDouble nz, RDouble &qx, RDouble &qy, RDouble &qz, RDouble *speciesFluxes)
{
    int nm = GetNSEquationNumber();
    int nTemperatureModel = GetTemperatureModel();
    int nSpecies = GetNumberOfSpecies();

    //! Allocate the memory space.
    RDouble *massFraction = new RDouble[nSpecies];
    RDouble *moleFraction = new RDouble[nSpecies];
    RDouble *speciesEnthalpy = new RDouble[nSpecies];
    //! It equals to the density multiplies the mass diffusion coefficient of the species, namely, rho * Ds.
    RDouble *speciesDiffusionCoef = new RDouble[nSpecies];
    RDouble *speciesVibrationEnergy = new RDouble[nSpecies];
    RDouble *speciesElectronEnergy = new RDouble[nSpecies];
    for (int s = 0; s < nSpecies; ++ s)
    {
        massFraction[s] = primitiveVars[nm + s];
    }

    //! Obtain the mole fraction of each species.
    ComputeMoleFractionByMassFraction(massFraction, moleFraction);

    //! Compute the mass diffusion coefficient of each species.
    ComputeSpeciesMassDiffusionCoefficient(massFraction, viscosityLaminar, viscosityTurbulence, speciesDiffusionCoef);

    //! Obtain the enthalpy of each species.
    GetEverySpeciesEnthalpy(transRotationTemperature, vibrationTemperatuare, electronTemperature, speciesEnthalpy);

    //! Compute the mass diffusion energies of the species, namely, rho * Ds * hs * dfsdn.
    RDouble diffusionCoef = 0.0;
    for (int s = 0; s < nSpecies; ++ s)
    {
        diffusionCoef = speciesDiffusionCoef[s] * speciesEnthalpy[s];
        qx += diffusionCoef * dfsdx[s];
        qy += diffusionCoef * dfsdy[s];
        qz += diffusionCoef * dfsdz[s];
    }

    RDouble vibrationEnergy[3] = {0.0, 0.0, 0.0};
    RDouble electronEnergy[3] = {0.0, 0.0, 0.0};
    RDouble peDivideDensity = GetElectronPressure(primitiveVars, electronTemperature) / MAX(primitiveVars[0], SMALL);

    if (nTemperatureModel > 1)
    {
        ComputeVibrationEnergy(vibrationTemperatuare, speciesVibrationEnergy);
        ComputeElectronEnergy(electronTemperature, speciesElectronEnergy);
        for (int s = 0; s < nSpecies; ++ s)
        {
            diffusionCoef = speciesDiffusionCoef[s] * speciesVibrationEnergy[s];
            vibrationEnergy[0] += diffusionCoef * dfsdx[s];
            vibrationEnergy[1] += diffusionCoef * dfsdy[s];
            vibrationEnergy[2] += diffusionCoef * dfsdz[s];

            diffusionCoef = speciesDiffusionCoef[s] * (speciesElectronEnergy[s] + peDivideDensity);
            electronEnergy[0] += diffusionCoef * dfsdx[s];
            electronEnergy[1] += diffusionCoef * dfsdy[s];
            electronEnergy[2] += diffusionCoef * dfsdz[s];
        }

        if (nTemperatureModel == 2)
        {
            speciesFluxes[nm + nSpecies] = nx * (vibrationEnergy[0] + electronEnergy[0]) + ny * (vibrationEnergy[1] + electronEnergy[1]) + nz * (vibrationEnergy[2] + electronEnergy[2]);
        }
        else
        {
            speciesFluxes[nm + nSpecies] = nx * vibrationEnergy[0] + ny * vibrationEnergy[1] + nz * vibrationEnergy[2];
            speciesFluxes[nm + nSpecies + 1] = nx * electronEnergy[0] + ny * electronEnergy[1] + nz * electronEnergy[2];
        }
    }

    //! Compute the mass diffusion flux of each species.
    for (int s = 0; s < nSpecies; ++ s)
    {
        speciesFluxes[nm + s] = speciesDiffusionCoef[s] * (nx * dfsdx[s] + ny * dfsdy[s] + nz * dfsdz[s]);
    }

    delete []massFraction;
    delete []moleFraction;
    delete []speciesEnthalpy;
    delete []speciesDiffusionCoef;
    delete []speciesVibrationEnergy;
    delete []speciesElectronEnergy;
}

void MixedGas::ComputeMixtureGasHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature,
                                       RDouble &transRotationConductivity, RDouble &vibrationConductivity, RDouble &electronConductivity)
{
    //int nTemperatureModel = GetTemperatureModel();
    if (ntmodel > 1)
    {
        //transRotationConductivity = GetMixedGasTranslationAndRotationHeatConductivity(primitiveVariables, transRotationTemperature, electronTemperature);
        //vibrationConductivity = GetMixedGasVibrationHeatConductivity(primitiveVariables, transRotationTemperature, vibrationTemperature, electronTemperature);
        //electronConductivity = GetMixedGasElectronHeatConductivity(primitiveVariables, transRotationTemperature, electronTemperature);
        ComputeMixedGasTranslationVibrationElectronHeatConductivity(primitiveVariables, transRotationTemperature, vibrationTemperature, electronTemperature, transRotationConductivity, vibrationConductivity, electronConductivity);
    }
    else     //! One-temperature model.
    {
        //! Compute heat conductivity of mixture gas.
        transRotationConductivity = ComputeMixtureGasHeatConductivityByWassilewaFormula(primitiveVariables, transRotationTemperature);
        vibrationConductivity = 0.0;
        electronConductivity = 0.0;
    }
}

void MixedGas::ComputeMixtureGasSpecificHeatAtConstantVolume(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature,
                                                   RDouble &cvtr, RDouble &cvv, RDouble &cve)
{
    RDouble speciesMassFraction[MAX_SPECIES_NUM];
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        speciesMassFraction[s] = primitiveVariables[nm + s];
    }
    if (ntmodel > 1)
    {
        cvtr = this->GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(speciesMassFraction);
        cvv  = this->GetMixedGasVibrationSpecificHeatAtConstantVolume(vibrationTemperature, speciesMassFraction);
        cve  = this->GetMixedGasElectronSpecificHeatAtConstantVolume(electronTemperature, speciesMassFraction);
    }
    else
    {
        this->ComputeConstantVolumeSpecificHeatByMassFraction(speciesMassFraction, transRotationTemperature, cvtr);
        cvv = 0.0;
        cve = 0.0;
    }
}

void MixedGas::ComputeEnergySourceDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                    RDouble ***chemicalSourceDerivatives, RDouble ***energySourceDerivatives, int nCellNumber)
{
    int nTemperatureModel = GetTemperatureModel();
    if (nTemperatureModel == 2)
    {
        ComputeEnergySourceDerivatives1(primitiveVariables, temperatures, cellVolume, chemicalSourceTerms,
            chemicalSourceDerivatives, energySourceDerivatives, nCellNumber);
    }
    else if (nTemperatureModel == 3)
    {
        //ComputeEnergySourceDerivatives1(primitiveVariables, temperatures, cellVolume, chemicalSourceTerms,
        //                                      chemicalSourceDerivatives, energySourceDerivatives, nCellNumber);
    }
}


void MixedGas::ComputeEnergySourceDiagonalDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                            RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber)
{
    //! Compute the partial derivative of energy source to the vibration and electron energy in the main diagonal position.
    int nTemperatureModel = GetTemperatureModel();
    if (nTemperatureModel == 2)
    {
        ComputeEnergySourceDiagonalDerivatives1(primitiveVariables, temperatures, cellVolume, chemicalSourceTerms, chemicalSourceDerivatives,
            energySourceDerivatives, nCellNumber);
    }
    else if (nTemperatureModel == 3)
    {
        ComputeEnergySourceDiagonalDerivatives2(primitiveVariables, temperatures, cellVolume, chemicalSourceTerms, chemicalSourceDerivatives,
            energySourceDerivatives, nCellNumber);
    }
}

RDouble MixedGas::GetMaximumSpeciesMassDiffusionCoef(RDouble *primitiveVars, RDouble viscosityLaminar, RDouble viscosityTurbulence)
{
    //! The function returns the maximum mass diffusion coefficient called (1 - Cs)/(1- Xs) * (mul/Scl + mut/Sct).
    RDouble maxDiffusionCoef = 0.0;
    if (numberOfSpecies < 1)
    {
        return maxDiffusionCoef;
    }

    //! Compute the mass diffusion coefficient of each species.
    GetSpeciesMassDiffusionCoef(primitiveVars, viscosityLaminar, viscosityTurbulence, workSpecies);

    //! Obtain the maximum value of Ds.
    for (int i = 0; i < nLastSpeciesIndex; ++ i)
    {
        maxDiffusionCoef = MAX(maxDiffusionCoef, workSpecies[i]);
    }

    return maxDiffusionCoef;     //! Returned value.
}

void MixedGas::GetSpeciesMassDiffusionCoef(RDouble *primitiveVars, RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble *speciesDiffusionCoef)
{
    //! The function returns the diffusion coefficients of each species which saved in rDs = (1 - Cs)/(1- Xs) * (mul/Scl + mut/Sct).
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVars[nm + i];
    }

    ComputeSpeciesMassDiffusionCoefficient(massFractions, viscosityLaminar, viscosityTurbulence, speciesDiffusionCoef);
}

void MixedGas::GetTemperatureAndCP(RDouble *prim,RDouble &Tm, RDouble &cp)
{
    int nchem = GetChemicalType();
    int nm = GetNSEquationNumber();

    RDouble gama = 1.4;     //! Specific heat ratio.
    RDouble &ro = prim[0];
    RDouble &pm = prim[4];
    RDouble R = GetCoefficientOfStateEquation();     //! The non-dimensional universal gas constant.
    cp = gama/(gama - 1.0) * R;    //! The constant pressure specific heat of perfect gas.
    Tm = pm/(ro * R);     //! Compute the temperature using the state equation of perfect gas.

    if (nchem == 1)     //! Mixture gas condition.
    {
        int ns = GetNumberOfSpecies();    //! Save the number of species.
        RDouble *fs = new RDouble[ns];
        for (int i = 0; i < ns; i++)
        {
            fs[i] = prim[i + nm];     //! Obtain the mass fraction of each species.
        }

        //! Obtain the average molecular weight of mixture gas.
        RDouble M = GetMixedGasMolecularWeight(fs);
        R = GetUniversalGasConstant();     //! The non-dimensional universal gas constant.
        //! Compute the temperature using the state equation of mixture gas.
        Tm = pm * M/(ro * R);
        
        //! Compute constant pressure specific heat of each species.
        ComputeConstantPressureSpecificHeatByMassFraction(fs, Tm, cp);

        delete []fs;     //! Remove the memory space.
    }
}

RDouble MixedGas::GetTemperature(RDouble *prim)
{
    int nchem = GetChemicalType();
    int nm = GetNSEquationNumber();
    //! Obtain the density and pressure.
    RDouble &ro = prim[0];
    RDouble &pm = prim[4];
    RDouble R = GetCoefficientOfStateEquation();
    //! Compute the temperature using the state equation of mixture gas.
    RDouble Tm = pm/(ro * R);

    if (nchem == 1)
    {
        //! Compute the temperature of mixture gas.
        int ns = GetNumberOfSpecies();     //! Save the number of species.
        RDouble *fs = new RDouble[ns];
        for (int i = 0; i < ns; i++)
        {
            fs[i] = prim[i + nm];     //! Obtain the mass fraction of each species.
        }

        //! Obtain the average molecular weight of mixture gas.
        RDouble M = GetMixedGasMolecularWeight(fs);
        R = GetUniversalGasConstant();     //! The non-dimensional universal gas constant.
        //! Compute the temperature using the state equation of mixture gas.
        Tm = pm * M/(ro * R);

        delete []fs;     //! Remove the memory space.
    }
    return Tm;
}

void MixedGas::GetEverySpeciesEnthalpy(RDouble trTemperature, RDouble vTemperature, RDouble eTemperature, RDouble *hs)
{
    if (ntmodel > 1)      //! Multi-temperature model.
    {
        ComputeSpeciesEnthalpy(trTemperature, vTemperature, eTemperature, hs);
    }
    else      //! Single temperature model.
    {
        ComputeSpeciesEnthalpy(trTemperature, hs);
    }
}

RDouble MixedGas::ComputeMolecularWeightReciprocalWithoutElectron(RDouble *prim, RDouble &ce_div_me)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for (int i = 0; i < numberOfSpecies; i++)
    {
        massFraction[i] = prim[i + nm];
    }

    RDouble M = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFraction, ce_div_me);
    return M;
}

RDouble MixedGas::GetMixedGasTranslationAndRotationEnergy(RDouble *prim, RDouble Ttr)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for(int i = 0; i < numberOfSpecies; i++)
    {
        massFraction[i] = prim[nm + i];
    }

    RDouble Etr = GetMixedGasTranslationAndRotationEnergy(Ttr, massFraction);
    return Etr;
}

RDouble MixedGas::GetMixedGasVibrationEnergy(RDouble *prim, RDouble Tv)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for(int i = 0; i < numberOfSpecies; ++ i)
    {
        massFraction[i] = prim[nm + i];
    }

    RDouble vibrationEnergy = GetMixedGasVibrationEnergy(Tv, massFraction);
    return vibrationEnergy;
}

RDouble MixedGas::GetMixedGasElectronEnergy(RDouble *prim, RDouble Te)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for(int i = 0; i < numberOfSpecies; ++ i)
    {
        massFraction[i] = prim[nm + i];
    }

    RDouble electronEnergy = GetMixedGasElectronEnergy(Te, massFraction);
    return electronEnergy;
}

RDouble MixedGas::GetMixedGasTransAndRotatSpecHeatAtConstVolume(RDouble *prim)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for(int i = 0; i < numberOfSpecies; i++)
    {
        massFraction[i] = prim[nm + i];
    }

    RDouble cvtr = GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(massFraction);
    return cvtr;
}

RDouble MixedGas::GetMixedGasVibrationSpecificHeatAtConstantVolume(RDouble* prim, RDouble Tv)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for(int i = 0; i < numberOfSpecies; i++)
    {
        massFraction[i] = prim[nm + i];
    }

    RDouble cvv = GetMixedGasVibrationSpecificHeatAtConstantVolume(Tv, massFraction);
    return cvv;
}

RDouble MixedGas::GetMixedGasElectronSpecificHeatAtConstantVolume(RDouble *prim, RDouble Te)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for(int i = 0; i < numberOfSpecies; i++)
    {
        massFraction[i] = prim[nm + i];
    }

    RDouble cve = GetMixedGasElectronSpecificHeatAtConstantVolume(Te, massFraction);
    return cve;
}

RDouble MixedGas::ComputeFrozenSoundSpeed(RDouble *primitiveVars)
{
    int nm = GetNSEquationNumber();
    int nChemical = GetChemicalType();
    int numberOfSpecies = GetNumberOfSpecies();
    RDouble massFraction[MAX_SPECIES_NUM];
    for(int s = 0; s < numberOfSpecies; ++ s)
    {
        massFraction[s] = primitiveVars[nm + s];
    }

    RDouble gasConstant = GetUniversalGasConstant();
    //! Obtain the specific heat at constant volume of translation and rotation mode.
    RDouble transRotationCv = GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(massFraction);

    //! Compute the coefficient namely alpha.
    RDouble massReciprocal = 0.0, ceDivideMe = 0.0;
    massReciprocal = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFraction, ceDivideMe);
    RDouble alpha = gasConstant * massReciprocal / transRotationCv;

    //! To compute the frozen sound speed.
    RDouble soundSpeed = sqrt((alpha + 1.0) * primitiveVars[IDX::IP] / primitiveVars[IDX::IR]);
    return soundSpeed;
}

void MixedGas::ComputeSpecificHeat(RDouble *primitiveVars, RDouble *temperatures, RDouble *spciesCvtr, RDouble *spciesCvv, RDouble &totalCvtr, RDouble &totalCvv)
{
    RDouble modifyCoef = 1.5;
    RDouble gasConstant = GetUniversalGasConstant();
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();

    //Compute the translational-rotational specific heat of species.
    this->ComputeTranslationAndRotationSpecificHeatAtConstantVolume(spciesCvtr);

    //Compute the vibrational specific heat of species.
    this->ComputeSpeciesConstantVolumeSpecificHeat(temperatures[IDX::ITV], spciesCvv);

    //Compute the total specific heat of mixture gas.
    totalCvtr = 0.0;
    totalCvv = 0.0;
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        totalCvtr += primitiveVars[nm + s] * spciesCvtr[s];

        spciesCvv[s] = modifyCoef * spciesCvv[s] - spciesCvtr[s] + (modifyCoef - 1.0) * gasConstant / speciesMass[s];
        totalCvv += primitiveVars[nm + s] * spciesCvv[s];
    }
}

void MixedGas::ComputeEnergy(RDouble *primitiveVars, RDouble *temperatures, RDouble *spciesEtr, RDouble *spciesEv, RDouble *spciesEnthalpy, RDouble &totalEtr, RDouble &totalEv, RDouble &totalEnthalpy)
{
    RDouble modifyCoef = 1.5;
    RDouble speciesCvtr[MAX_SPECIES_NUM];
    //Obtain the translation-rotation specific heat at constant volume.
    ComputeTranslationAndRotationSpecificHeatAtConstantVolume(speciesCvtr);

    //Obtain the thermodynamic properties of species.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    //Compute the translational-rotational energies of species.
    for (int s = 0; s < this->numberOfSpecies; ++ s)
    {
        spciesEtr[s]  = speciesCvtr[s] * temperatures[IDX::ITT];
        spciesEtr[s] += thermodynamicProperties[s].hs0 / referenceEnergy;
    }

    RDouble gasConstant = this->GetUniversalGasConstant();// dimensional universal gas constant.
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();

    //Compute the enthalpies of species.
    this->ComputeSpeciesEnthalpy(temperatures[IDX::ITV], spciesEnthalpy);

    //Compute the vibrational energies of species.
    totalEtr = 0.0;
    totalEv = 0.0;
    totalEnthalpy = 0.0;
    for (int s = 0; s < this->numberOfSpecies; ++ s)
    {
        spciesEv[s] = modifyCoef * spciesEnthalpy[s] - gasConstant * temperatures[IDX::ITV] / speciesMass[s];
        spciesEv[s] = spciesEv[s] - speciesCvtr[s] * temperatures[IDX::ITV];
        spciesEv[s] = spciesEv[s] - thermodynamicProperties[s].hs0 / referenceEnergy;

        spciesEnthalpy[s] = spciesEtr[s] + spciesEv[s] + gasConstant * temperatures[IDX::ITT] / speciesMass[s];

        totalEtr += primitiveVars[nm + s] * spciesEtr[s];
        totalEv  += primitiveVars[nm + s] * spciesEv [s];
        totalEnthalpy += primitiveVars[nm + s] * spciesEnthalpy[s];
    }
}

void MixedGas::ComputePressureDerivatives(RDouble *primitiveVars, RDouble *temperatures, RDouble &alpha, RDouble &beta, RDouble &eta, RDouble *speciesBeta)
{
    using namespace IDX;
    RDouble sigma1 = 0.0, sigma2 = 0.0;
    //Obtain the thermodynamic properties of species.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    int *speciesCharge = this->molecularProperty->GetIonTypeOfSpecies();
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();
    RDouble tmp = speciesCharge[nLastSpeciesIndex] * this->electronMass / speciesMass[nLastSpeciesIndex];
    RDouble squareVelocity = primitiveVars[IU] * primitiveVars[IU] + primitiveVars[IV] * primitiveVars[IV] + primitiveVars[IW] * primitiveVars[IW];
    sigma1 = 1.0 / (1.0 + tmp);
    sigma2 = tmp * sigma1;

    RDouble electronMassFraction = 0.0, electronEv = 0.0;
    if (this->nElectronIndex >= 0)
    {
        electronMassFraction = speciesMass[this->nElectronIndex];
    }

    RDouble gasConstant = this->GetUniversalGasConstant();
    RDouble totalRMass = this->ComputeMolecularWeightReciprocal(primitiveVars);
    RDouble speciesCvtr[MAX_SPECIES_NUM], speciesCvv[MAX_SPECIES_NUM], totalCvtr, totalCvv;
    this->ComputeSpecificHeat(primitiveVars, temperatures, speciesCvtr, speciesCvv, totalCvtr, totalCvv);
    alpha = (gasConstant / totalCvtr) * (totalRMass - electronMassFraction / this->electronMass);
    eta   = (gasConstant / totalCvv)  * (electronMassFraction / this->electronMass);

    RDouble speciesEtr[MAX_SPECIES_NUM], speciesEv[MAX_SPECIES_NUM], speciesHs[MAX_SPECIES_NUM];
    RDouble totalEtr, totalEv, totalH;
    this->ComputeEnergy(primitiveVars, temperatures, speciesEtr, speciesEv, speciesHs, totalEtr, totalEv, totalH);
    if (this->nElectronIndex >= 0)
    {
        electronEv = speciesEv[this->nElectronIndex];
    }
    tmp = sigma1 * temperatures[ITT] / speciesMass[nLastSpeciesIndex] + sigma2 * temperatures[ITV] / this->electronMass;
    beta = alpha * (0.5 * squareVelocity - sigma1 * speciesEtr[nLastSpeciesIndex]) + gasConstant * tmp;
    beta -= eta * (sigma1 * speciesEv[nLastSpeciesIndex] + sigma2 * electronEv);

    for (int j = 0; j < nLastSpeciesIndex; ++ j)
    {
        tmp = speciesCharge[j] * this->electronMass / speciesMass[j];
        speciesBeta[j]  = alpha * (sigma1 * (1.0 + tmp) * speciesEtr[nLastSpeciesIndex] - speciesEtr[j]);
        speciesBeta[j] += eta * (sigma1 * speciesEv[nLastSpeciesIndex] + sigma2 * electronEv - speciesEv[j]);
        speciesBeta[j] += eta * sigma1 * tmp * (speciesEv[nLastSpeciesIndex] - electronEv);
        speciesBeta[j] += gasConstant * temperatures[ITT] * (1.0 / speciesMass[j] - sigma1 * (1.0 + tmp) / speciesMass[nLastSpeciesIndex]);
        speciesBeta[j] += (sigma1 * tmp - sigma2) * gasConstant * temperatures[ITV] / this->electronMass;
    }
    speciesBeta[nLastSpeciesIndex] = 0.0;
    if (this->nElectronIndex >= 0)
    {
        speciesBeta[nElectronIndex] = 0.0;
    }
}

void MixedGas::ModifyVibrationEnergy(char *filepath)
{
    FILE *pfile = fopen(filepath, "wb");

    RDouble Tv = 0.0, Ttr = 0.0;
    RDouble deltaTv = 100.0, tmp1 = 0.0, eps = 0.0, delta_eps = 0.0, hiTv = 0.0, cpTv = 0.0, cvvTv = 0.0;

    RDouble speciesCvtr[MAX_SPECIES_NUM], spciesEtr[MAX_SPECIES_NUM], spciesEnthalpy[MAX_SPECIES_NUM], spciesEv[MAX_SPECIES_NUM];
    RDouble speciesEv1[MAX_SPECIES_NUM], speciesEe[MAX_SPECIES_NUM], speciesCvv[MAX_SPECIES_NUM], speciesCve[MAX_SPECIES_NUM];
    //Obtain the thermodynamic properties of species.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    //Obtain the translation-rotation specific heat at constant volume.
    this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCvtr);

    RDouble gasConstant = rjmk;// dimensional universal gas constant.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    RDouble refTemperature = referenceParameterDimensional.GetTemperature();

    while (Tv <= 30000.0)
    {
        Ttr = Tv;
        //Compute the translational-rotational energies of species.
        for (int s = 0; s < this->numberOfSpecies; ++ s)
        {
            spciesEtr[s] = speciesCvtr[s] * Ttr + thermodynamicProperties[s].hs0;
        }

        //Compute the enthalpies of species.
        this->ComputeSpeciesEnthalpyDimensional(Tv / refTemperature, spciesEnthalpy);

        this->ComputeDimensionalVibrationEnergy(Tv / refTemperature, speciesEv1);
        this->ComputeDimensionalElectronEnergy(Tv / refTemperature, speciesEe);
        this->ComputeDimensionalVibrationSpecificHeatAtConstantVolume(Tv / refTemperature, speciesCvv);
        this->ComputeDimensionalElectronSpecificHeatAtConstantVolume(Tv / refTemperature, speciesCve);

        //Compute the vibrational energies of species.
        fprintf(pfile, "%f", Tv);
        for (int s = 0; s < this->numberOfSpecies; ++ s)
        {
            tmp1 = gasConstant * Tv / speciesMass[s] + speciesCvtr[s] * Tv + thermodynamicProperties[s].hs0;
            spciesEv[s] = spciesEnthalpy[s] - tmp1;
            speciesEv1[s] += speciesEe[s];

            delta_eps = speciesEv1[s] - spciesEv[s];
            hiTv = speciesEv1[s] + tmp1;
            eps = hiTv / spciesEnthalpy[s];

            cvvTv = speciesCvv[s] + speciesCve[s];
            cpTv = cvvTv + gasConstant / speciesMass[s] + speciesCvtr[s];
            //fprintf(pfile, " %.6e %.6e %.6e", delta_eps, eps, hiTv);
            //fprintf(pfile, " %.8e %.8e %.8e", cvvTv, cpTv, hiTv);
            fprintf(pfile, " %.8e", cvvTv);
            //fprintf(pfile, " %.8e", speciesEv1[s]);
        }
        fprintf(pfile, "\r\n");

        Tv += deltaTv;
    }


    fclose(pfile);
}

void MixedGas::ComputeSpecificHeatByFitMethod(RDouble temperature, RDouble *cvSpecies, bool isVibrationModel/* = true*/)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    RDouble referenceTemperature = referenceParameterDimensional.GetTemperature();
    RDouble *coef, cp;
    RDouble T1, T2, T3, T4;
    T1 = temperature * referenceTemperature;
    //T1 = MAX(T1, 50.0);
    T2 = T1 * T1;
    T3 = T1 * T2;
    T4 = T1 * T3;

    if (isVibrationModel)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType < 2)
            {
                cvSpecies[ispecies] = 0.0;
                continue;
            }
            coef = thermodynamicManager->GetEnthalpyFitCoefficient(ispecies);
            cp = coef[1] + 2 * coef[2] * T1 + 3 * coef[3] * T2 + 4 * coef[4] * T3 + 5 * coef[5] * T4;
            cvSpecies[ispecies] = (cp - 1.0) * gasConstant / speciesMass[ispecies];
        }
    }
    else
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType >= 2)
            {
                cvSpecies[ispecies] = 0.0;
                continue;
            }
            coef = thermodynamicManager->GetEnthalpyFitCoefficient(ispecies);
            cp = coef[1] + 2 * coef[2] * T1 + 3 * coef[3] * T2 + 4 * coef[4] * T3 + 5 * coef[5] * T4;
            cvSpecies[ispecies] = (cp - 1.0) * gasConstant / speciesMass[ispecies];
        }
    }
}

void MixedGas::ComputeEnthalpyByFitMethod(RDouble temperature, RDouble *hSpecies, bool isVibrationModel/* = true*/)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    RDouble referenceTemperature = referenceParameterDimensional.GetTemperature();
    RDouble *coef, hs;
    RDouble T1, T2, T3, T4, T5;
    T1 = temperature * referenceTemperature;
    //T1 = MAX(T1, 50.0);
    T2 = T1 * T1;
    T3 = T1 * T2;
    T4 = T1 * T3;
    T5 = T1 * T4;

    if (isVibrationModel)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType < 2)
            {
                hSpecies[ispecies] = 0.0;
                continue;
            }
            coef = thermodynamicManager->GetEnthalpyFitCoefficient(ispecies);
            hs = coef[0] + coef[1] * T1 + coef[2] * T2 + coef[3] * T3 + coef[4] * T4 + coef[5] * T5;
            hSpecies[ispecies] = hs * gasConstant / speciesMass[ispecies];
        }
    }
    else
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType >= 2)
            {
                hSpecies[ispecies] = 0.0;
                continue;
            }
            coef = thermodynamicManager->GetEnthalpyFitCoefficient(ispecies);
            hs = coef[0] + coef[1] * T1 + coef[2] * T2 + coef[3] * T3 + coef[4] * T4 + coef[5] * T5;
            hSpecies[ispecies] = hs * gasConstant / speciesMass[ispecies];
        }
    }
}

void MixedGas::ComputeSpecificHeatByPolynomialFitMethod(RDouble temperature, RDouble *cvSpecies, bool isVibrationModel/* = true*/)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    RDouble referenceTemperature = referenceParameterDimensional.GetTemperature();
    RDouble cp, T1, T2, T3, T4;
    T1 = temperature * referenceTemperature;
    //T1 = MAX(T1, 1.0);
    T2 = T1 * T1;
    T3 = T1 * T2;
    T4 = T1 * T3;

    //Obtain the interpolation coefficients for the polynomial.
    this->GetLinearlyAveragedPolynomialCoefficientsT2(T1, this->polynomialCoefT2);

    if (isVibrationModel)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType < 2)
            {
                cvSpecies[ispecies] = 0.0;
                continue;
            }
            cp = polynomialCoefT2[ispecies][1] + 2 * polynomialCoefT2[ispecies][2] * T1 + 3 * polynomialCoefT2[ispecies][3] * T2 + 4 * polynomialCoefT2[ispecies][4] * T3 + 5 * polynomialCoefT2[ispecies][5] * T4;
            cvSpecies[ispecies] = (cp - 1.0) * gasConstant / speciesMass[ispecies];
        }
    }
    else
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType >= 2)
            {
                cvSpecies[ispecies] = 0.0;
                continue;
            }
            cp = polynomialCoefT2[ispecies][1] + 2 * polynomialCoefT2[ispecies][2] * T1 + 3 * polynomialCoefT2[ispecies][3] * T2 + 4 * polynomialCoefT2[ispecies][4] * T3 + 5 * polynomialCoefT2[ispecies][5] * T4;
            cvSpecies[ispecies] = (cp - 1.0) * gasConstant / speciesMass[ispecies];
        }
    }
}

void MixedGas::ComputeEnthalpyByPolynomialFitMethod(RDouble temperature, RDouble *hSpecies, bool isVibrationModel/* = true*/)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    RDouble referenceTemperature = referenceParameterDimensional.GetTemperature();
    RDouble hs, T1, T2, T3, T4, T5;
    T1 = temperature * referenceTemperature;
    //T1 = MAX(T1, 1.0);
    T2 = T1 * T1;
    T3 = T1 * T2;
    T4 = T1 * T3;
    T5 = T1 * T4;

    //Obtain the interpolation coefficients for the polynomial.
    this->GetLinearlyAveragedPolynomialCoefficientsT2(T1, this->polynomialCoefT2);

    if (isVibrationModel)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType < 2)
            {
                hSpecies[ispecies] = 0.0;
                continue;
            }
            hs = polynomialCoefT2[ispecies][0] + polynomialCoefT2[ispecies][1] * T1 + polynomialCoefT2[ispecies][2] * T2 + polynomialCoefT2[ispecies][3] * T3 + polynomialCoefT2[ispecies][4] * T4 + polynomialCoefT2[ispecies][5] * T5;
            hSpecies[ispecies] = hs * gasConstant / speciesMass[ispecies];
        }
    }
    else
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType >= 2)
            {
                hSpecies[ispecies] = 0.0;
                continue;
            }
            hs = polynomialCoefT2[ispecies][0] + polynomialCoefT2[ispecies][1] * T1 + polynomialCoefT2[ispecies][2] * T2 + polynomialCoefT2[ispecies][3] * T3 + polynomialCoefT2[ispecies][4] * T4 + polynomialCoefT2[ispecies][5] * T5;
            hSpecies[ispecies] = hs * gasConstant / speciesMass[ispecies];
        }
    }
}

void ChemCompressData(DataContainer *&cdata)
{
    gas->CompressData(cdata);
}

void ChemDecompressData(DataContainer *cdata)
{
    gas->DecompressData(cdata);
}


void  computeOutPressure(RDouble &Pressure, RDouble &staticPressure,  RDouble &staticPressureRelaxCorf )
{   
    Pressure= MAX(  staticPressureRelaxCorf * staticPressure + ( 1.0 - staticPressureRelaxCorf )*Pressure ,  staticPressure  );

//  cout <<"  staticPressure=   "<<  staticPressure  <<"  staticPressureRelaxCorf=   "<<  staticPressureRelaxCorf  << endl;
//  abort();
}




void MixedGas::SetMaximumSpecies(RDouble *fs)
{
    RDouble *Ms = this->molecularProperty->GetMolecularWeightDimensional();
    string *species_name = this->molecularProperty->GetNameOfSpecies();

    char    elementName[4]   = {'O', 'N', 'C', 'H'}; 
    RDouble elementWeight[4] = {1.599940E-2, 1.400674E-2, 1.400674E-2, 1.0E-3};
    RDouble elementMaximum[4] = {0.0, 0.0, 0.0, 0.0};

    int coef = 0;
    RDouble elementMaximumC1;

    RDouble *xs = new RDouble[this->numberOfSpecies]();

    for (int i = 0; i < this->numberOfSpecies; i++)
    {
        maximumSpecies[i] = 1.0;    //ini
        xs[i] = fs[i] / Ms[i];
    }

    if(nSpeciesLimit==0)return;

    //compute the maximum of each element
    for (int n = 0; n < 4; n++)
    {
        for (int i = 0; i < this->numberOfSpecies; i++)
        {
            if (FindChemicalElement(species_name[i], elementName[n], coef))
            {
                elementMaximum[n] += xs[i] * coef * elementWeight[n];
            }
        }
    }

    //compute the maximum of each Specie
    for (int i = 0; i < this->numberOfSpecies; i++)
    {
        for (int n = 0; n < 4; n++)
        {
            if (FindChemicalElement(species_name[i], elementName[n], coef))
            {
                elementMaximumC1 =  Ms[i] * elementMaximum[n] / (elementWeight[n] * coef) ;
                maximumSpecies[i]=MIN( maximumSpecies[i] , elementMaximumC1 );
            }
        }
    }
    delete [] xs;    xs = NULL;


}

}
}
