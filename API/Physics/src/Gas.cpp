#include "Gas.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "TK_Exit.h"


namespace PHSPACE
{
namespace GAS_SPACE
{
    Gas *gas;

    int GetIndexBySpeciesName(string speciesName)
    {
        int nIndex = -1;
        if (speciesName == "O")
        {
            nIndex = 0;
        }
        else if (speciesName == "O2")
        {
            nIndex = 1;
        }
        else if (speciesName == "NO")
        {
            nIndex = 2;
        }
        else if (speciesName == "N")
        {
            nIndex = 3;
        }
        else if (speciesName == "N2")
        {
            nIndex = 4;
        }
        else if (speciesName == "O+")
        {
            nIndex = 5;
        }
        else if (speciesName == "O2+")
        {
            nIndex = 6;
        }
        else if (speciesName == "NO+")
        {
            nIndex = 7;
        }
        else if (speciesName == "N+")
        {
            nIndex = 8;
        }
        else if (speciesName == "N2+")
        {
            nIndex = 9;
        }
        else if (speciesName == "e-")
        {
            nIndex = 10;
        }
        else if (speciesName == "C")
        {
            nIndex = 11;
        }
        else if (speciesName == "C2")
        {
            nIndex = 12;
        }
        else if (speciesName == "CO")
        {
            nIndex = 13;
        }
        else if (speciesName == "CO2")
        {
            nIndex = 14;
        }
        else if (speciesName == "CN")
        {
            nIndex = 15;
        }
        else
        {
            nIndex = -1;
        }
        return nIndex;
    }

    string GetSpeciesNameByIndex(int nIndex)
    {
        string speciesName = "";
        if (nIndex == 0)
        {
            speciesName = "O";
        }
        else if (nIndex == 1)
        {
            speciesName = "O2";
        }
        else if (nIndex == 2)
        {
            speciesName = "NO";
        }
        else if (nIndex == 3)
        {
            speciesName = "N";
        }
        else if (nIndex == 4)
        {
            speciesName = "N2";
        }
        else if (nIndex == 5)
        {
            speciesName = "O+";
        }
        else if (nIndex == 6)
        {
            speciesName = "O2+";
        }
        else if (nIndex == 7)
        {
            speciesName = "NO+";
        }
        else if (nIndex == 8)
        {
            speciesName = "N+";
        }
        else if (nIndex == 9)
        {
            speciesName = "N2+";
        }
        else if (nIndex == 10)
        {
            speciesName = "e-";
        }
        else if (nIndex == 11)
        {
            speciesName = "C";
        }
        else if (nIndex == 12)
        {
            speciesName = "C2";
        }
        else if (nIndex == 13)
        {
            speciesName = "CO";
        }
        else if (nIndex == 14)
        {
            speciesName = "CO2";
        }
        else if (nIndex == 15)
        {
            speciesName = "CN";
        }
        else
        {
            speciesName = "";
        }
        return speciesName;
    }

    void SpeciesNameToInteger(int numberOfSpecies, string *nameList, Int1D *speciesOrder)
    {
        //Initialization.
        for (int n = 0; n < numberOfSpecies; ++ n)
        {
            (*speciesOrder)(n) = -1;
        }
        for (int n = 0; n < numberOfSpecies; ++ n)
        {
            int nIndex = GetIndexBySpeciesName(nameList[n]);
            if (nIndex >= 0)
            {
                (*speciesOrder)(nIndex) = n;
            }
        }
    }
    void SpeciesNameToInteger(int numberOfSpecies, string *nameList, int *speciesOrder)
    {
        //Initialization.
        for (int n = 0; n < numberOfSpecies; ++ n)
        {
            speciesOrder[n] = -1;
        }
        for (int n = 0; n < numberOfSpecies; ++ n)
        {
            int nIndex = GetIndexBySpeciesName(nameList[n]);
            if (nIndex >= 0)
            {
                speciesOrder[nIndex] = n;
            }
        }
    }

    bool IsEqualGasModel(Int1D *speciesOrder1, Int1D *speciesOrder2)
    {
        //the index of the array are O, O2, NO, N, N2, O+, O2+, NO+, N+, N2+, e-, C, C2, CO, CO2, CN, respectively.
        //the elements indicate the order of the species in the gas model, -1 indicates the species does not exist.
        int totalNumber = 16;
        for (int n = 0; n < totalNumber; ++ n)
        {
            if ((*speciesOrder1)(n) != (*speciesOrder2)(n))
            {
                return false;
            }
            //if (((*speciesOrder1)(n) >= 0 && (*speciesOrder2)(n) < 0) ||
            //    ((*speciesOrder1)(n) < 0 && (*speciesOrder2)(n) >= 0))
            //{
            //    return false;
            //}
        }
        return true;
    }

    Gas::Gas()
    {
        wallMultiTemperature = 0;
        if (GlobalDataBase::IsExist("wallMultiTemperature", PHINT, 1))
        {
            wallMultiTemperature = GlobalDataBase::GetIntParaFromDB("wallMultiTemperature");
        }

        if (GlobalDataBase::IsExist("useHyflowSetting", PHINT, 1))
        {
            if( GlobalDataBase::GetIntParaFromDB("useHyflowSetting") == 1)
            {
                wallMultiTemperature = 1;
            }
        }
    };

#ifdef USE_HYFLOW_CODE
void Gas::GetAirInfo_Modify(const RDouble &h1, RDouble &t, RDouble &p, RDouble &den, RDouble &a)
{
    //! USA 1976 Model for atmosphere.
    using namespace GAS_SPACE;

    RDouble refMolecularWeight = GlobalDataBase::GetDoubleParaFromDB("refMolecularWeight");
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble r = 1000.0 * GeneralGasConstant / refMolecularWeight;      //be careful

    RDouble h    = h1 * 1000.0;
    RDouble g0   = 9.80665;
    RDouble rp   = 6.356766e6;
    RDouble g    = (rp / (rp + h)) * (rp / (rp + h)) * g0;
    RDouble t0   = 288.15;
    RDouble p0   = 101325.0;
    RDouble ato1 = -0.0065;
    RDouble ato2 = 0.001;

    RDouble rho = 0.0;

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


    if (h <= 11019.0)
    {
        RDouble al1 = (t11 - t0)/11019.0;
        t   = t0 + al1 * h;
        p   = p0 * pow(t/t0, -g/(r * al1));
        p   = p * 100.0;
        //rho = rho0 * pow(t/t0, -1.0 - g/(r*al1));
    }
    else if (h <= 20063.0)
    {
        t   = t11;
        p   = p11   * exp(-g*(h-11019.0)/(r*t11));
        p   = p * 100.0;
        //rho = rho11 * exp(-g*(h-11019.0)/(r*t11));
    }
    else if (h <= 32162.0)
    {
        RDouble al2 = (t32-t20)/(32162.0-20063.0);
        t   = t11+al2*(h-20063.0);
        p   = p20 * pow(t/t11, -g/(r*al2));
        p   = p * 100.0;
        //rho = rho20 * pow(t/t11, -1.0-g/(r*al2));
    }
    else if (h <= 47350.0)
    {
        RDouble al3 = (t47-t32)/(47350.0-32162.0);
        t   = t32+al3*(h-32162.0);
        p   = p32 * pow(t/t32, -g/(r*al3));
        p   = p * 100.0;
        //rho = rho32 * pow(t/t32, -1.0-g/(r*al3));
    }
    else if (h <= 51000.0)
    {
        t   = t47;
        p   = p47 * exp(-g*(h-47350.0)/(r*t47));
        p   = p * 100.0;
        //rho = rho47*exp(-g*(h-47350.0)/(r*t47));
    }
    else if(h > 51000 && h <= 61000)
    {
        t = 270.65 - 0.0026 * (h-51000);
        p = 70.73853 * pow((1 - 1.034546 * 1e-5 * (h - 51000)), 12.20115);
        //rho = 9.116011 * 1e-4 * pow((1 - 1.034546 * 1e-5 * (h - 51000)), 11.20115);
    }
    else if(h > 61000 && h <= 71000)
    {
        t = 270.65 - 0.00265 * (h-51000);
        p = 73.73853 * pow((1 - 1.034546 * 1e-5 * (h - 51000)), 12.20115);
        //rho = 9.416011 * 1e-4 * pow((1 - 1.034546 * 1e-5 * (h - 51000)), 11.20115);
    }
    else if(h > 71000 && h < 86000)
    {
        t = 216.65 - 0.00199 * (h - 71000);
        p = 4.658392 * pow((1 - 9.317494 * 1e-6 * (h - 71000)), 17.08161);
        //rho = 7.521057 * 1e-5 * pow((1 - 9.317494 * 1e-6 * (h - 71000)), 16.08161);
    }
    else if(h >= 86000 && h <= 91000)
    {
        t = 186.87;
        p = 1.61325 * 0.1 * (2.273 + 1.042 * 0.000001 * h) * exp((86000 - h) / 5470);
        //rho = 7.0833475 * 0.000001 * exp((86000 - h) / 5470);
    }
    else if(h > 91000 && h <= 100000)
    {
        t = 1082.36 - 0.01961 * h + 1.07366e-7 * h * h;
        p = 10.81744-0.000212186 * h + 0.00000000104341 * h * h;
    }
    else
    {
        TK_Exit::ExceptionExit("Error: height must be less than 100km!!!\n", true);
    }

    rho = p / (r * t);
    a   = sqrt(refGama * r * t);
    den = rho;
}
#endif

}

}