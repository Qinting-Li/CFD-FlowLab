#include "PHMpi.h"
#include "LESSolverStruct.h"
#include "GradientOperation.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "Param_LESSolverStruct.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"

using namespace std;

namespace PHSPACE
{

LESSolverStruct::LESSolverStruct()
{

}

LESSolverStruct::~LESSolverStruct()
{
    DeAllocateGlobalVariables();
    FreeControlParameters();
}

void LESSolverStruct::AllocateGlobalVar(Grid *gridIn)
{
    //LESSolver::InitParameter();
    Param_LESSolverStruct * parameters = GetControlParameters();

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range IFACE, JFACE, KFACE;
    GetRange(ni, nj, nk, 0, 0, IFACE, JFACE, KFACE);

    Range Direction(1, 3);
    Range N(1, 6);
    Range M(1, 3);

    RDouble5D *subgridStress = new RDouble5D(IFACE, JFACE, KFACE, Range(1, 6), Direction, fortranArray);
    *subgridStress = 0.0;
    grid->UpdateDataPtr("subgridStress", subgridStress);

    RDouble5D *subgridHeatFlux = new RDouble5D(IFACE, JFACE, KFACE, Range(1, 3), Direction, fortranArray);
    *subgridHeatFlux = 0.0;
    grid->UpdateDataPtr("subgridHeatFlux", subgridHeatFlux);

    RDouble3D *cellLengthScale = new RDouble3D(I, J, K, fortranArray);
    *cellLengthScale = 0.0;
    grid->UpdateDataPtr("cellLengthScale", cellLengthScale);
    
    RDouble4D *faceLengthScale = new RDouble4D(IFACE, JFACE, KFACE, Direction ,fortranArray);
    *faceLengthScale = 0.0;
    grid->UpdateDataPtr("faceLengthScale", faceLengthScale);

    //RDouble4D *wallFunction = new RDouble4D(IFACE, JFACE, KFACE, Direction ,fortranArray);
    RDouble3D *wallFunction = new RDouble3D(I, J, K, fortranArray);
    *wallFunction = 1.0;
    grid->UpdateDataPtr("wallFunction", wallFunction);

    //RDouble4D *subgridScaleEnergy = new RDouble4D(I, J, K, Direction, fortranArray);
    //*subgridScaleEnergy = 0.0;
    //grid->UpdateDataPtr("kSGS", subgridScaleEnergy);

    RDouble3D *anisotropicConstant = new RDouble3D(I, J, K, fortranArray);
    *anisotropicConstant = 0.0;
    grid->UpdateDataPtr("anisotropicConstant", anisotropicConstant);

    RDouble3D *isotropicConstant = new RDouble3D(I, J, K, fortranArray);
    *isotropicConstant = 0.0;
    grid->UpdateDataPtr("isotropicConstant", isotropicConstant);

    RDouble3D *turbulentPrandtlNumber = new RDouble3D(I, J, K, fortranArray);
    *turbulentPrandtlNumber = 0.9;
    grid->UpdateDataPtr("turbulentPrandtlNumber", turbulentPrandtlNumber);

    RDouble4D *strainRateTensor = new RDouble4D(I, J, K, Range(0, 6), fortranArray);
    *strainRateTensor = 0.0;
    grid->UpdateDataPtr("strainRateTensor", strainRateTensor);

    //RDouble3D *strainRateMagnitude = new RDouble3D(I, J, K, fortranArray);
    //*strainRateMagnitude = 0.0;
    //grid->UpdateDataPtr("strainRateMagnitude", strainRateMagnitude);

    int subgridScaleModel = parameters->GetSubgridScaleModel();
    if (subgridScaleModel == 2 || subgridScaleModel == 4)
    {
        RDouble4D *leonardStress = new RDouble4D(I, J, K, N, fortranArray);
        RDouble4D *modeledStress = new RDouble4D(I, J, K, N, fortranArray);
        RDouble4D *leonardTemperatureFlux = new RDouble4D(I, J, K, Range(0, 2), fortranArray);
        RDouble4D *modeledTemperatureFlux = new RDouble4D(I, J, K, Range(0, 2), fortranArray);
        RDouble4D *testFilteredQ = new RDouble4D(I, J, K, Range(0, 3), fortranArray);
        //RDouble3D *testFilteredDensity = new RDouble3D(I, J, K, fortranArray);
        RDouble3D *testFilteredTemperature = new RDouble3D(I, J, K, fortranArray);
        RDouble4D *testFilteredStrainRateTensor = new RDouble4D(I, J, K, Range(0, 6), fortranArray);
        RDouble4D *testFilteredRhoUiUj = new RDouble4D(I, J, K, N, fortranArray);
        RDouble4D *testFilteredRhoUT = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *testFilteredAlphaIJ = new RDouble4D(I, J, K, Range(0, 6), fortranArray);
        RDouble4D *betaIJ = new RDouble4D(I, J, K, N, fortranArray);
        RDouble4D *gradT = new RDouble4D(I, J, K, Range(0, 2), fortranArray);
        RDouble4D *testFilteredGradT = new RDouble4D(I, J, K, Range(0, 2), fortranArray);
        RDouble4D *testFilteredRhoStrainRateMagnitudeGradT = new RDouble4D(I, J, K, Range(0, 2), fortranArray);

        *leonardStress = 0.0;
        *modeledStress= 0.0;
        *leonardTemperatureFlux = 0.0;
        *modeledTemperatureFlux = 0.0;
        *testFilteredQ = 0.0;
        //*testFilteredDensity = 0.0;
        *testFilteredTemperature = 0.0;
        *testFilteredStrainRateTensor = 0.0;
        *testFilteredRhoUiUj = 0.0;
        *testFilteredRhoUT = 0.0;
        *testFilteredAlphaIJ = 0.0;
        *betaIJ = 0.0;
        *gradT = 0.0;
        *testFilteredGradT = 0.0;
        *testFilteredRhoStrainRateMagnitudeGradT = 0.0;
 
        grid->UpdateDataPtr("leonardStress", leonardStress);
        grid->UpdateDataPtr("modeledStress", modeledStress);
        grid->UpdateDataPtr("leonardTemperatureFlux", leonardTemperatureFlux);
        grid->UpdateDataPtr("modeledTemperatureFlux", modeledTemperatureFlux);
        grid->UpdateDataPtr("testFilteredQ", testFilteredQ);
        //grid->UpdateDataPtr("testFilteredDensity", testFilteredDensity);
        grid->UpdateDataPtr("testFilteredTemperature", testFilteredTemperature);
        grid->UpdateDataPtr("testFilteredStrainRateTensor", testFilteredStrainRateTensor);
        grid->UpdateDataPtr("testFilteredRhoUiUj", testFilteredRhoUiUj);
        grid->UpdateDataPtr("testFilteredRhoUT", testFilteredRhoUT);
        grid->UpdateDataPtr("testFilteredAlphaIJ", testFilteredAlphaIJ);
        grid->UpdateDataPtr("betaIJ", betaIJ);
        grid->UpdateDataPtr("gradT", gradT);
        grid->UpdateDataPtr("testFilteredGradT", testFilteredGradT);
        grid->UpdateDataPtr("testFilteredRhoStrainRateMagnitudeGradT", testFilteredRhoStrainRateMagnitudeGradT);
    }

    Init(grid);
}

void LESSolverStruct::DeAllocateGlobalVar(Grid *gridIn)
{
    Param_LESSolverStruct * parameters = GetControlParameters();
    StructGrid *grid = StructGridCast(gridIn);

    RDouble5D *subgridStress = reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridStress"));
    RDouble5D *subgridHeatFlux = reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridHeatFlux"));
    delete subgridStress;
    delete subgridHeatFlux;

    RDouble3D *cellLengthScale = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cellLengthScale"));
    RDouble4D *faceLengthScale = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("faceLengthScale"));
    RDouble3D *wallFunction = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));
    delete cellLengthScale;
    delete faceLengthScale;
    delete wallFunction;

    //RDouble4D *subgridScaleEnergy = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("kSGS"));
    //delete subgridScaleEnergy;

    RDouble3D *anisotropicConstant = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D *isotropicConstant = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));
    RDouble3D *turbulentPrandtlNumber = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    RDouble4D *strainRateTensor = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));

    delete isotropicConstant;
    delete anisotropicConstant;
    delete turbulentPrandtlNumber;
    delete strainRateTensor;

    int subgridScaleModel = parameters->GetSubgridScaleModel();
    if (subgridScaleModel == 2 || subgridScaleModel == 4)
    { 
        RDouble4D *leonardStress = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardStress"));
        RDouble4D *modeledStress = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledStress"));
        RDouble4D *leonardTemperatureFlux = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardTemperatureFlux"));
        RDouble4D *modeledTemperatureFlux = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledTemperatureFlux"));
        RDouble4D *testFilteredQ = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));
        //RDouble3D *testFilteredDensity = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredDensity"));
        RDouble3D *testFilteredTemperature = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredTemperature"));
        RDouble4D *testFilteredStrainRateTensor = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredStrainRateTensor"));
        RDouble4D *testFilteredRhoUiUj = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUiUj"));
        RDouble4D *testFilteredRhoUT = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUT"));
        RDouble4D *testFilteredAlphaIJ = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredAlphaIJ"));
        RDouble4D *betaIJ = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("betaIJ"));
        RDouble4D *gradT = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradT"));
        RDouble4D *testFilteredGradT = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredGradT"));
        RDouble4D *testFilteredRhoStrainRateMagnitudeGradT = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoStrainRateMagnitudeGradT"));

        delete leonardStress;
        delete modeledStress;
        delete leonardTemperatureFlux;
        delete modeledTemperatureFlux;
        delete testFilteredQ;
        //delete testFilteredDensity;
        delete testFilteredTemperature;
        delete testFilteredStrainRateTensor;
        delete testFilteredRhoUiUj;
        delete testFilteredRhoUT;
        delete testFilteredAlphaIJ;
        delete betaIJ;
        delete gradT;
        delete testFilteredGradT;
        delete testFilteredRhoStrainRateMagnitudeGradT;
    }
}

LIB_EXPORT void LESSolverStruct::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_LESSolverStruct();
    controlParameters->Init();
}

LIB_EXPORT Param_LESSolverStruct * LESSolverStruct::GetControlParameters()
{
    return static_cast<Param_LESSolverStruct *> (controlParameters);
}

void LESSolverStruct::Init(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ifvfd = GlobalDataBase::GetIntParaFromDB("ifvfd");
    if(ifvfd == 0)
    {
        ComputeCellLengthScale(grid);
    }
    else
    {
        ComputeCellLengthScaleforStructHighOrder(grid);
    }
    //ComputeFaceLengthScale(grid);
}

void LESSolverStruct::ComputeCellLengthScale(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D &cellLengthScale = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("cellLengthScale"));

    RDouble4D &area = *(grid->GetFaceArea());
    RDouble3D &vol  = *(grid->GetCellVolume());

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int nDim = GetDim();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    //double areaL[3], areaR[3], delta[3];

    int i, j, k;

    //! inner cells
    for (k = kst; k <= ked; ++ k)
    {
        for (j = jst; j <= jed; ++ j)
        {
            for (i = ist; i <= ied; ++ i)
            {
                RDouble len3 = 1.0;

                for (int iSurface = 1; iSurface <= nDim; ++iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    //areaL[iSurface - 1] = area(i, j, k, iSurface);
                    //areaR[iSurface - 1] = area(il, jl, kl, iSurface);
                    RDouble areaL = area(i, j, k, iSurface);
                    RDouble areaR = area(il, jl, kl, iSurface);

                    RDouble volume = vol(i, j, k);

                    //delta[iSurface - 1] = 2.0 * volume / (areaL[iSurface - 1] + areaR[iSurface - 1]);
                    RDouble delta = 2.0 * volume / (areaL + areaR);

                    len3 = len3 * delta;
                }
                
                cellLengthScale(i, j, k) = pow(len3, third);
            }
        }
    }

    //! ghost cells, for i = 0 and ni;
    i = 0;

    for (k = kst; k <= ked; ++ k)
    {
        for (j = jst; j <= jed; ++ j)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    i = ni;

    for (k = kst; k <= ked; ++ k)
    {
        for (j = jst; j <= jed; ++ j)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    //! ghost cells, for j = 0 and nj;
    j = 0;

    for (k = kst; k <= ked; ++ k)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    j = nj;

    for (k = kst; k <= ked; ++ k)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    //! ghost cells, for k = 0 and nk;
    k = 0;

    for (j = jst; j <= jed; ++ j)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    k = nk;

    for (j = jst; j <= jed; ++ j)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }
}

void LESSolverStruct::ComputeCellLengthScaleforStructHighOrder(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D &cellLengthScale = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("cellLengthScale"));

    RDouble4D &area = *(grid->GetFaceArea());
    RDouble3D &vol  = *(grid->GetCellVolume());

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int nDim = GetDim();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    //double areaL[3], areaR[3], delta[3];

    int i, j, k;

    //! inner cells
    for (k = kst; k <= ked; ++ k)
    {
        for (j = jst; j <= jed; ++ j)
        {
            for (i = ist; i <= ied; ++ i)
            {
                RDouble len3 = 1.0;

                for (int iSurface = 1; iSurface <= nDim; ++iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    //areaL[iSurface - 1] = area(i, j, k, iSurface);
                    //areaR[iSurface - 1] = area(il, jl, kl, iSurface);
                    RDouble areaL = area(i, j, k, iSurface);
                    RDouble areaR = area(il, jl, kl, iSurface);

                    RDouble volume = vol(i, j, k);

                    //delta[iSurface - 1] = 2.0 * volume / (areaL[iSurface - 1] + areaR[iSurface - 1]);
                    RDouble delta = 2.0 * volume / (areaL + areaR);

                    len3 = len3 * delta;
                }
                
                cellLengthScale(i, j, k) = pow(len3, third);
            }
        }
    }

    GhostCell3D(cellLengthScale, ni, nj, nk);

    RDouble a = cellLengthScale(-1, 1, 1);
    /*
    //! ghost cells, for i = 0 and ni;
    i = 0;

    for (k = kst; k <= ked; ++ k)
    {
        for (j = jst; j <= jed; ++ j)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = 0 + il1;
                int jl = j + jl1;
                int kl = k + kl1;
                int il2 = 0 + il1 * 2;
                int jl2 = j + jl1 * 2;
                int kl2 = k + kl1 * 2;

                RDouble areaL = area(0-il1, j-jl1, k-kl1, iSurface);
                RDouble areaM = area(0    , j    , k    , iSurface);
                RDouble areaR = area(0+il1, j+jl1, k+kl1, iSurface);

                RDouble volumeL = vol(-1, j, k);
                RDouble volumeR = vol(0 , j, k);

                RDouble deltaL = 2.0 * volumeL / (areaL + areaM);

                len = len * deltaL;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    i = ni;

    for (k = kst; k <= ked; ++ k)
    {
        for (j = jst; j <= jed; ++ j)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    //! ghost cells, for j = 0 and nj;
    j = 0;

    for (k = kst; k <= ked; ++ k)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    j = nj;

    for (k = kst; k <= ked; ++ k)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    //! ghost cells, for k = 0 and nk;
    k = 0;

    for (j = jst; j <= jed; ++ j)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    k = nk;

    for (j = jst; j <= jed; ++ j)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }
    */
}

void LESSolverStruct::ComputeFaceLengthScale(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &faceLengthScale = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("faceLengthScale"));
    RDouble3D &cellLengthScale = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cellLengthScale"));

    int nDim = GetDim();

    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int ist, ied, jst, jed, kst, ked;
        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble cellLengthL = cellLengthScale(i, j, k);
                    RDouble cellLengthR = cellLengthScale(il, jl, kl);

                    RDouble faceLength = 0.5 * (cellLengthScale(i, j, k) + cellLengthScale(il, jl, kl));

                    faceLengthScale(il, jl, kl, iSurface) = faceLength;
                }
            }
        }
    }

}

void LESSolverStruct::ComputeWallFunction(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    
    Param_LESSolverStruct * parameters = GetControlParameters();
    int wallDampingFunctionType = parameters->GetWallDampingFunctionType();

    RDouble3D &wallFunction = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    if (wallDampingFunctionType == 0)
    {
        wallFunction = 1.0;
    }
    else if (wallDampingFunctionType == 1)
    {
        WallFunctionofVanDriest(grid);
    }
    else if (wallDampingFunctionType == 2)
    {
        WallFunctionofDXB(grid);
    }
    else if (wallDampingFunctionType == 3)
    {
        WallFunctionofPiomelli(grid);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("wallDampingFunctionType", wallDampingFunctionType);
    }
}

/*
void LESSolverStruct::WallFunctionofVanDriestOld(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D & wallDistance = * StructGridCast(grid)->GetWallDist();

    RDouble4D & q = * reinterpret_cast<RDouble4D *> ( grid->GetDataPtr("q") );

    RDouble3D & visl = * reinterpret_cast<RDouble3D *> ( grid->GetDataPtr("visl") );

    RDouble5D &gradPrimtiveVarFaceXofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceXofLES"));
    RDouble5D &gradPrimtiveVarFaceYofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceYofLES"));
    RDouble5D &gradPrimtiveVarFaceZofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceZofLES"));

    RDouble4D &wallFunction = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("wallFunction"));

    Param_LESSolverStruct * parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();

    int nDim = GetDim();

    const double kappa = 0.41;
    const double Ap = 25.0;

    //! i direction;
    int iSurface = 1;

    int ist, ied, jst, jed, kst, ked;
    grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

    int il1, jl1, kl1;
    GetNsurfIndex(iSurface, il1, jl1, kl1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                int indexI = static_cast<int>(abs(i-ist-0.5) + ist + 0.5);
                RDouble wallDistL = wallDistance(indexI, j, k);

                int indexIl = static_cast<int>(ied + 0.5 - abs(il-ied-0.5));
                RDouble wallDistR = wallDistance(indexIl, jl, kl);

                RDouble wallDist = 0.5 * (wallDistL + wallDistR);

                RDouble mul = 0.5 * (visl(i, j, k) + visl(il, jl, kl));

                RDouble rho = 0.5 * (q(i, j, k, IR) + q(il, jl, kl, IR));

                RDouble dudx = gradPrimtiveVarFaceXofLES(i, j, k, IU, iSurface);
                RDouble dudy = gradPrimtiveVarFaceYofLES(i, j, k, IU, iSurface);
                RDouble dudz = gradPrimtiveVarFaceZofLES(i, j, k, IU, iSurface);

                RDouble dvdx = gradPrimtiveVarFaceXofLES(i, j, k, IV, iSurface);
                RDouble dvdy = gradPrimtiveVarFaceYofLES(i, j, k, IV, iSurface);
                RDouble dvdz = gradPrimtiveVarFaceZofLES(i, j, k, IV, iSurface);

                RDouble dwdx = gradPrimtiveVarFaceXofLES(i, j, k, IW, iSurface);
                RDouble dwdy = gradPrimtiveVarFaceYofLES(i, j, k, IW, iSurface);
                RDouble dwdz = gradPrimtiveVarFaceZofLES(i, j, k, IW, iSurface);

                RDouble s11 = dudx;
                RDouble s22 = dvdy;
                RDouble s33 = dwdz;
                RDouble s12 = half * (dudy + dvdx);
                RDouble s13 = half * (dudz + dwdx);
                RDouble s23 = half * (dvdz + dwdy);

                RDouble sijsij = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                RDouble ss = sqrt(two * sijsij);

                RDouble yplus = ss * (kappa * wallDist) * (kappa * wallDist) * rho / mul * refReNumber;

                wallFunction(il, jl, kl, iSurface) = 1.0 - exp( -yplus / Ap );
            }
        }
    }

    //! j direction
    iSurface = 2;

    grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

    GetNsurfIndex(iSurface, il1, jl1, kl1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                int indexJ = static_cast<int>(abs(j-jst-0.5) + jst + 0.5);
                RDouble wallDistL = wallDistance(i, indexJ, k);

                int indexJl = static_cast<int>(jed + 0.5 - abs(jl-jed-0.5));
                RDouble wallDistR = wallDistance(il, indexJl, kl);

                RDouble wallDist = 0.5 * (wallDistL + wallDistR);

                RDouble mul = 0.5 * (visl(i, j, k) + visl(il, jl, kl));

                RDouble rho = 0.5 * (q(i, j, k, IR) + q(il, jl, kl, IR));

                RDouble dudx = gradPrimtiveVarFaceXofLES(i, j, k, IU, iSurface);
                RDouble dudy = gradPrimtiveVarFaceYofLES(i, j, k, IU, iSurface);
                RDouble dudz = gradPrimtiveVarFaceZofLES(i, j, k, IU, iSurface);

                RDouble dvdx = gradPrimtiveVarFaceXofLES(i, j, k, IV, iSurface);
                RDouble dvdy = gradPrimtiveVarFaceYofLES(i, j, k, IV, iSurface);
                RDouble dvdz = gradPrimtiveVarFaceZofLES(i, j, k, IV, iSurface);

                RDouble dwdx = gradPrimtiveVarFaceXofLES(i, j, k, IW, iSurface);
                RDouble dwdy = gradPrimtiveVarFaceYofLES(i, j, k, IW, iSurface);
                RDouble dwdz = gradPrimtiveVarFaceZofLES(i, j, k, IW, iSurface);

                RDouble s11 = dudx;
                RDouble s22 = dvdy;
                RDouble s33 = dwdz;
                RDouble s12 = half * (dudy + dvdx);
                RDouble s13 = half * (dudz + dwdx);
                RDouble s23 = half * (dvdz + dwdy);

                RDouble sijsij = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                RDouble ss = sqrt(two * sijsij);

                RDouble yplus = ss * (kappa * wallDist) * (kappa * wallDist) * rho / mul * refReNumber;

                wallFunction(il, jl, kl, iSurface) = 1.0 - exp(-yplus / Ap);
            }
        }
    }

    //! k direction
    iSurface = 3;

    grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

    GetNsurfIndex(iSurface, il1, jl1, kl1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                int indexK = static_cast<int>(abs(k-kst-0.5) + kst + 0.5);
                RDouble wallDistL = wallDistance(i, j, indexK);

                int indexKl = static_cast<int>(ked + 0.5 - abs(kl-ked-0.5));
                RDouble wallDistR = wallDistance(il, jl, indexKl);

                RDouble wallDist = 0.5 * (wallDistL + wallDistR);

                RDouble mul = 0.5 * (visl(i, j, k) + visl(il, jl, kl));

                RDouble rho = 0.5 * (q(i, j, k, IR) + q(il, jl, kl, IR));

                RDouble dudx = gradPrimtiveVarFaceXofLES(i, j, k, IU, iSurface);
                RDouble dudy = gradPrimtiveVarFaceYofLES(i, j, k, IU, iSurface);
                RDouble dudz = gradPrimtiveVarFaceZofLES(i, j, k, IU, iSurface);

                RDouble dvdx = gradPrimtiveVarFaceXofLES(i, j, k, IV, iSurface);
                RDouble dvdy = gradPrimtiveVarFaceYofLES(i, j, k, IV, iSurface);
                RDouble dvdz = gradPrimtiveVarFaceZofLES(i, j, k, IV, iSurface);

                RDouble dwdx = gradPrimtiveVarFaceXofLES(i, j, k, IW, iSurface);
                RDouble dwdy = gradPrimtiveVarFaceYofLES(i, j, k, IW, iSurface);
                RDouble dwdz = gradPrimtiveVarFaceZofLES(i, j, k, IW, iSurface);

                RDouble s11 = dudx;
                RDouble s22 = dvdy;
                RDouble s33 = dwdz;
                RDouble s12 = half * (dudy + dvdx);
                RDouble s13 = half * (dudz + dwdx);
                RDouble s23 = half * (dvdz + dwdy);

                RDouble sijsij = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                RDouble ss = sqrt(two * sijsij);

                RDouble yplus = ss * (kappa * wallDist) * (kappa * wallDist) * rho / mul * refReNumber;

                wallFunction(il, jl, kl, iSurface) = 1.0 - exp(-yplus / Ap);
            }
        }
    }
}
*/

void LESSolverStruct::WallFunctionofVanDriest(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D & wallDistance = * StructGridCast(grid)->GetWallDist();

    RDouble4D & q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D & visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D &wallFunction = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));
    //RDouble3D &strainRateMagnitude = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("strainRateMagnitude"));

    Param_LESSolverStruct * parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();

    const double kappa = 0.41;
    const double Ap = 25.0;

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble wallDist = wallDistance(i, j, k);

                RDouble mul = visl(i, j, k);

                RDouble rho = q(i, j, k, IR);

                RDouble s11 = strainRateTensor(i, j, k, 1);
                RDouble s22 = strainRateTensor(i, j, k, 2);
                RDouble s33 = strainRateTensor(i, j, k, 3);
                RDouble s12 = strainRateTensor(i, j, k, 4);
                RDouble s13 = strainRateTensor(i, j, k, 5);
                RDouble s23 = strainRateTensor(i, j, k, 6);

                RDouble ss = strainRateTensor(i, j, k, 0);

                RDouble yplus = ss * (kappa * wallDist) * (kappa * wallDist) * rho / mul * refReNumber;

                wallFunction(i, j, k) = 1.0 - exp(-yplus / Ap);
            }
        }
    }

    GhostCell3D(wallFunction, ni, nj, nk);
}

void LESSolverStruct::WallFunctionofDXB(Grid *gridIn)
{
    using namespace IDX;

    StructGrid * grid = StructGridCast(gridIn);

    RDouble3D & wallDistance = * StructGridCast(grid)->GetWallDist();

    RDouble4D & q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D & visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));
    //RDouble3D &strainRateMagnitude = * reinterpret_cast<RDouble3D *> ( grid->GetDataPtr("strainRateMagnitude") );

    Param_LESSolverStruct * parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble rgam = parameters->GetRgam();

    const double kappa = 0.41;
    const double Ap = 25.0;
    const double yploge = 1000.0;

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble wallDist = wallDistance(i, j, k);

                RDouble mul = visl(i, j, k);

                RDouble rho = q(i, j, k, IR);

                RDouble s11 = strainRateTensor(i, j, k, 1);
                RDouble s22 = strainRateTensor(i, j, k, 2);
                RDouble s33 = strainRateTensor(i, j, k, 3);
                RDouble s12 = strainRateTensor(i, j, k, 4);
                RDouble s13 = strainRateTensor(i, j, k, 5);
                RDouble s23 = strainRateTensor(i, j, k, 6);

                RDouble ss = strainRateTensor(i, j, k, 0);

                RDouble yplus = ss * (kappa * wallDist) * (kappa * wallDist) * rho / mul * refReNumber;

                RDouble ffd = tanh((yploge/yplus) * (yploge/yplus) * (yploge/yplus));
                RDouble vanDriestInner = rgam * MIN((1.0 - exp(-pow(yplus/Ap, 7.0/6.0))), (1.0 - exp(-yplus / Ap)));
                RDouble vanDriestOut = 1.0;

                wallFunction(i, j, k) = vanDriestInner * ffd + vanDriestOut * (1.0 - ffd);
            }
        }
    }

    GhostCell3D(wallFunction, ni, nj, nk);
}

void LESSolverStruct::WallFunctionofPiomelli(Grid *gridIn)
{
    using namespace IDX;

    StructGrid * grid = StructGridCast(gridIn);

    RDouble3D & wallDistance = * StructGridCast(grid)->GetWallDist();

    RDouble4D & q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D & visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));
    //RDouble3D &strainRateMagnitude = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("strainRateMagnitude"));

    Param_LESSolverStruct * parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();

    const double kappa = 0.41;
    const double Ap = 25.0;

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble wallDist = wallDistance(i, j, k);

                RDouble mul = visl(i, j, k);

                RDouble rho = q(i, j, k, IR);

                RDouble s11 = strainRateTensor(i, j, k, 1);
                RDouble s22 = strainRateTensor(i, j, k, 2);
                RDouble s33 = strainRateTensor(i, j, k, 3);
                RDouble s12 = strainRateTensor(i, j, k, 4);
                RDouble s13 = strainRateTensor(i, j, k, 5);
                RDouble s23 = strainRateTensor(i, j, k, 6);

                RDouble ss = strainRateTensor(i, j, k, 0);

                RDouble yplus = ss * (kappa * wallDist) * (kappa * wallDist) * rho / mul * refReNumber;

                wallFunction(i, j, k) = sqrt(1.0 - exp(-(yplus/Ap) * (yplus/Ap) * (yplus/Ap)));
            }
        }
    }

    GhostCell3D(wallFunction, ni, nj, nk);
}

void LESSolverStruct::ComputeViscousCoeff(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    ComputeStrainRateTensor(grid);

    ComputeWallFunction(grid);

    Param_LESSolverStruct * parameters = GetControlParameters();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    if (subgridScaleModel == 1)
    {
        Smagorinsky(grid);
    }
    else if (subgridScaleModel == 2)
    {
        DynamicSmagViscosityCompressible(grid);
    }
    else if (subgridScaleModel == 3)
    {
        WALE(grid);
    }
    else if (subgridScaleModel == 4)
    {
        DynamicSmagViscosityIncompressible(grid);
    }
    else if (subgridScaleModel == 5)
    {
        DynamicSmagViscosity(grid);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }
}

void LESSolverStruct::SmagorinskyOld(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &faceLengthScale = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("faceLengthScale"));
    RDouble4D &wallFunction = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble5D &subgridStress = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridStress"));
    //RDouble4D &subgridEnergy = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("kSGS"));
    RDouble5D &subgridHeatFlux = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridHeatFlux"));

    RDouble5D &gradPrimtiveVarFaceXofLES = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceXofLES"));
    RDouble5D &gradPrimtiveVarFaceYofLES = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceYofLES"));
    RDouble5D &gradPrimtiveVarFaceZofLES = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceZofLES"));

    RDouble5D &gradtx = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceXofLES"));
    RDouble5D &gradty = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceYofLES"));
    RDouble5D &gradtz = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceZofLES"));

    Param_LESSolverStruct * parameters = GetControlParameters();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    RDouble refReNumber = parameters->GetRefReNumber();

    RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

    RDouble smagConstant = parameters->GetSmagConstant();
    RDouble isotropicConstant = parameters->GetIsotropicConstant();

    int turbViscousCutType = parameters->GetTurbViscousCutType();

    RDouble refGama = parameters->GetRefGama();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);

    double vistmax  = 0.0;
    double vistmin  = 1.0e30;
    double logsmall = log(SMALL);

    //int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    //ComputeSubgridScaleEnergy(grid);

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    int imin = 1;
    int jmin = 1;
    int kmin = 1;

    int nDim = GetDim();

    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int ist, ied, jst, jed, kst, ked;
        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble deltaBar = faceLengthScale(il ,jl ,kl, iSurface);

                    RDouble dudx = gradPrimtiveVarFaceXofLES(il, jl, kl, IU, iSurface);
                    RDouble dudy = gradPrimtiveVarFaceYofLES(il, jl, kl, IU, iSurface);
                    RDouble dudz = gradPrimtiveVarFaceZofLES(il, jl, kl, IU, iSurface);

                    RDouble dvdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IV, iSurface);

                    RDouble dwdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IW, iSurface);

                    RDouble dtdx = gradtx(il, jl, kl, ITT, iSurface);
                    RDouble dtdy = gradty(il, jl, kl, ITT, iSurface);
                    RDouble dtdz = gradtz(il, jl, kl, ITT, iSurface);

                    RDouble s11 = dudx;
                    RDouble s22 = dvdy;
                    RDouble s33 = dwdz;
                    RDouble s12 = half * (dudy + dvdx);
                    RDouble s13 = half * (dudz + dwdx);
                    RDouble s23 = half * (dvdz + dwdy);

                    RDouble sijsij = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                    RDouble ss = sqrt(two * sijsij);
                    RDouble skk = s11 + s22 + s33;

                    RDouble rho = 0.5 * (q(i, j, k, IR) + q(il, jl, kl, IR));
                    RDouble mul = 0.5 * (visl(i, j, k) + visl(il, jl, kl));
                    RDouble dampingFunction = wallFunction(il, jl, kl, iSurface);

                    RDouble mut = (smagConstant * deltaBar * dampingFunction) * (smagConstant * deltaBar * dampingFunction) * rho * ss * refReNumber;

                    if (turbViscousCutType == 1)
                    {
                        mut = MAX(mut - mul, 0.0);
                    }
                    else if (turbViscousCutType == 2)
                    {
                        mut = MAX(mut, 0.0);
                    }
                    else if (turbViscousCutType == 3)
                    {
                        mut = MAX(mut, -mul);
                    }

                    mut = MIN(eddyViscosityLimit, mut);

                    RDouble kSGS = 2.0 * isotropicConstant * (deltaBar * dampingFunction) * (deltaBar * dampingFunction) * rho * ss * ss * refReNumber;

                    subgridStress(il, jl, kl, 1, iSurface) = -2.0 * mut * (s11 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    subgridStress(il, jl, kl, 2, iSurface) = -2.0 * mut * (s22 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    subgridStress(il, jl, kl, 3, iSurface) = -2.0 * mut * (s33 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    subgridStress(il, jl, kl, 4, iSurface) = -2.0 * mut * s12;
                    subgridStress(il, jl, kl, 5, iSurface) = -2.0 * mut * s13;
                    subgridStress(il, jl, kl, 6, iSurface) = -2.0 * mut * s23;

                    subgridHeatFlux(il, jl, kl, 1, iSurface) = -mut * specificHeatAtConstantPressure * oPrandtlTurbulence * dtdx;
                    subgridHeatFlux(il, jl, kl, 2, iSurface) = -mut * specificHeatAtConstantPressure * oPrandtlTurbulence * dtdy;
                    subgridHeatFlux(il, jl, kl, 3, iSurface) = -mut * specificHeatAtConstantPressure * oPrandtlTurbulence * dtdz;

                    if (vistmax < mut)
                    {
                        vistmax = mut;
                        imax = i;
                        jmax = j;
                        kmax = k;
                    }

                    if (vistmin > mut)
                    {
                        vistmin = mut;
                        imin = i;
                        jmin = j;
                        kmin = k;
                    }

                }
            }
        }
    }

    double vist_max = vistmax;
    double vist_min = vistmin;

    grid->UpdateData("vist_max", &vist_max, PHDOUBLE, 1);
    grid->UpdateData("vist_min", &vist_min, PHDOUBLE, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorVistmax)
        {
            cout << "vistmax = " << vistmax << " " << imax << " " << jmax << " " << kmax << "\n";
            cout << "vistmin = " << vistmin << " " << imin << " " << jmin << " " << kmin << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void LESSolverStruct::Smagorinsky(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D &cellLengthScale = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cellLengthScale"));
    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));

    Param_LESSolverStruct * parameters = GetControlParameters();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    RDouble refReNumber  = parameters->GetRefReNumber();

    RDouble smagConstant = parameters->GetSmagConstant();

    RDouble iConstant = parameters->GetIsotropicConstant();

    double vistmax  = 0.0;
    double vistmin  = 1.0e30;
    double logsmall = log(SMALL);

    //int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    //ComputeSubgridScaleEnergy(grid);
    int ifvfd = GlobalDataBase::GetIntParaFromDB("ifvfd");

    int ist, ied, jst, jed, kst, ked;
    if(ifvfd == 0)
    {
        grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    }
    else
    {
        grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);
    }
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble deltaBar = cellLengthScale(i, j, k);

                RDouble dampingFunction = wallFunction(i, j, k);

                //RDouble mut = (smagConstant * deltaBar * dampingFunction) * (smagConstant * deltaBar * dampingFunction) * rho * ss * refReNumber;

                RDouble cr = (smagConstant * deltaBar * dampingFunction) * (smagConstant * deltaBar * dampingFunction);

                RDouble ci = iConstant * (deltaBar * dampingFunction) * (deltaBar * dampingFunction);

                anisotropicConstant(i, j, k) = cr;

                isotropicConstant(i, j, k) = ci;
            }
        }
    }
    
    RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();

    turbulentPrandtlNumber = prandtlTurbulence;
}

/*
void LESSolverStruct::SmagorinskyViscosity(Grid *grid_in)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(grid_in);

    RDouble4D &xfn  = *(grid->GetFaceNormalX());
    RDouble4D &yfn  = *(grid->GetFaceNormalY());
    RDouble4D &zfn  = *(grid->GetFaceNormalZ());
    RDouble4D &area = *(grid->GetFaceArea());
    RDouble3D &vol  = *(grid->GetCellVolume());

    RDouble4D & q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D & visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble5D & subgridStress = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridStress"));

    RDouble4D & gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D & gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D & gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    RDouble3D & wallDistance = * StructGridCast(grid)->GetWallDist();

    Param_LESSolverStruct * parameters = GetControlParameters();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    RDouble refReNumber  = parameters->GetRefReNumber();

    int wallDampingFunctionType = parameters->GetWallDampingFunctionType();

    RDouble yploge = 0.0;
    if (wallDampingFunctionType == 2)
    {
        yploge = parameters->GetYploge();
    }

    RDouble smagConstant = parameters->GetSmagConstant();

    int deltaFunctionType = parameters->GetDeltaFunctionType();

    RDouble rgam = parameters->GetRgam();

    int turbViscousCutType = parameters->GetTurbViscousCutType();

    double vistmax  = 0.0;
    double vistmin  = 1.0e30;
    double logsmall = log(SMALL);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble areai = 0.5 * (area(i, j, k, 1) + area(i+1, j, k, 1));
                RDouble areaj = 0.5 * (area(i, j, k, 2) + area(i, j+1, k, 2));
                RDouble areak = 0.5 * (area(i, j, k, 3) + area(i, j, k+1, 3));
                RDouble deltai = vol(i, j, k) / areai;
                RDouble deltaj = vol(i, j, k) / areaj;
                RDouble deltak = vol(i, j, k) / areak;

                RDouble deltaBar = 0.0;;
                GetDeltaBar(deltaFunctionType, deltai, deltaj, deltak, deltaBar);

                RDouble dudx = gradUVWTCellCenterX(i, j, k, 0);
                RDouble dudy = gradUVWTCellCenterY(i, j, k, 0);
                RDouble dudz = gradUVWTCellCenterZ(i, j, k, 0);

                RDouble dvdx = gradUVWTCellCenterX(i, j, k, 1);
                RDouble dvdy = gradUVWTCellCenterY(i, j, k, 1);
                RDouble dvdz = gradUVWTCellCenterZ(i, j, k, 1);

                RDouble dwdx = gradUVWTCellCenterX(i, j, k, 2);
                RDouble dwdy = gradUVWTCellCenterY(i, j, k, 2);
                RDouble dwdz = gradUVWTCellCenterZ(i, j, k, 2);

                RDouble s11 = dudx;
                RDouble s22 = dvdy;
                RDouble s33 = dwdz;
                RDouble s12 = half * (dudy + dvdx);
                RDouble s13 = half * (dudz + dwdx);
                RDouble s23 = half * (dvdz + dwdy);

                RDouble sij2 = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                RDouble ss = sqrt(two * sij2);
                RDouble rho = q(i, j, k, IR);
                RDouble mul = visl(i, j, k);

                //! Calculate damping function fVD.
                RDouble Ap = 25.0;
                RDouble rKappa = 0.41;
                RDouble wallDist = wallDistance(i, j, k);
                RDouble rKappa_walldist = rKappa * wallDist;
                RDouble rKappa_walldist2 = rKappa_walldist * rKappa_walldist;
                RDouble yp = ss * (rKappa * wallDist) * (rKappa * wallDist) * rho / mul * refReNumber;
                RDouble fVD = 1.0;
                if (wallDampingFunctionType == 0)
                {
                    fVD = 1.0;
                }
                else if (wallDampingFunctionType == 1)
                {
                    fVD = 1.0 - exp(-yp / Ap);
                }
                else if (wallDampingFunctionType == 2)
                {
                    //! Developed by Deng XiaoBing.
                    RDouble ffd = tanh((yploge/yp) * (yploge/yp) * (yploge/yp));
                    RDouble rVDinner = rgam * MIN((1.0 - exp(-pow(yp/Ap, 7.0/6.0))), (1.0 - exp(-yp / Ap)));
                    RDouble rVDout = 1.0;
                    fVD = rVDinner * ffd + rVDout * (1.0 - ffd);
                }
                else if (wallDampingFunctionType == 3)
                {
                    //! Piomelli 1993.
                    fVD = sqrt(1.0 - exp(-(yp/Ap) * (yp/Ap) * (yp/Ap)));
                }
                else
                {
                    TK_Exit::UnexpectedVarValue("wallDampingFunctionType", wallDampingFunctionType);
                }

                //! Calculate mut.
                RDouble csmag_fVD_delta = smagConstant * fVD * deltaBar;
                RDouble csmag_fVD_delta2 = csmag_fVD_delta * csmag_fVD_delta;
                RDouble mut = rho * csmag_fVD_delta2 * ss * refReNumber;

                if (turbViscousCutType == 1)
                {
                    mut = MAX(mut - mul, 0.0);
                }
                else if (turbViscousCutType == 2)
                {
                    mut = MAX(mut, 0.0);
                }

                mut = MIN(eddyViscosityLimit, mut);

                if (vistmax < vist(i,j,k))
                {
                    vistmax = vist(i,j,k);
                    imax = i;
                    jmax = j;
                    kmax = k;
                }
            }
        }
    }

    double vist_max = vistmax;
    double vist_min = vistmin;

    grid->UpdateData("vist_max", &vist_max, PHDOUBLE, 1);
    grid->UpdateData("vist_min", &vist_min, PHDOUBLE, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorVistmax)
        {
            cout << "vistmax = " << vistmax << " " << imax << " " << jmax << " " << kmax << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif

}
*/

void LESSolverStruct::DynamicSmagViscosity(Grid *grid_in)
{
    using namespace IDX;

    RDouble u1, u2, u3;
    RDouble s11, s22, s33, s12, s13, s23;
    RDouble sij2, smod;
    RDouble dudx, dudy, dudz;
    RDouble dvdx, dvdy, dvdz;
    RDouble dwdx, dwdy, dwdz;
    RDouble rho, mut, mul;

    RDouble rL11, rL22, rL33, rL12, rL13, rL23;
    RDouble rM11, rM22, rM33, rM12, rM13, rM23;
    RDouble tfss, tfsijtfsij;

    StructGrid *grid = StructGridCast(grid_in);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();
    int nis1 = ni - 1;
    int njs1 = nj - 1;
    int nks1 = nk - 1;

    RFloat4D & q    = * reinterpret_cast<RFloat4D *> (grid->GetDataPtr("q"  ));
    RFloat3D & visl = * reinterpret_cast<RFloat3D *> (grid->GetDataPtr("visl"));
    RFloat3D & vist = * reinterpret_cast<RFloat3D *> (grid->GetDataPtr("vist"));

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble refReNumber = parameters->GetRefReNumber();

    RDouble rfwr = parameters->GetTestFilterScale();

    int turbViscousCutType = parameters->GetTurbViscousCutType();

    int *filterDirection = parameters->GetFilterDirection();

    int *averageDirection = parameters->GetAverageDirection();

    RDouble vistmax  = 0.0;
    RDouble vistmin  = 1.0e30;

    Range I(1, ni-1);
    Range J(1, nj-1);
    Range K(1, nk-1);
    if (nk == 1) K.setRange(1, 1);

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    Range II(-1, ni+1);
    Range JJ(-1, nj+1);
    Range KK(-1, nk+1);

    RFloat3D gradux(II, JJ, KK, fortranArray);
    RFloat3D graduy(II, JJ, KK, fortranArray);
    RFloat3D graduz(II, JJ, KK, fortranArray);

    RFloat3D gradvx(II, JJ, KK, fortranArray);
    RFloat3D gradvy(II, JJ, KK, fortranArray);
    RFloat3D gradvz(II, JJ, KK, fortranArray);

    RFloat3D gradwx(II, JJ, KK, fortranArray);
    RFloat3D gradwy(II, JJ, KK, fortranArray);
    RFloat3D gradwz(II, JJ, KK, fortranArray);

    GradCenter(grid, q, gradux, graduy, graduz, IU);
    GradCenter(grid, q, gradvx, gradvy, gradvz, IV);
    GradCenter(grid, q, gradwx, gradwy, gradwz, IW);

    Range III(0, ni);
    Range JJJ(0, nj);
    Range KKK(0, nk);

    RFloat3D *tfu1, *tfu2, *tfu3;
    RFloat3D *tfuifuj11, *tfuifuj22, *tfuifuj33, *tfuifuj12, *tfuifuj13, *tfuifuj23;
    RFloat3D *tfs11, *tfs22, *tfs33, *tfs12, *tfs13, *tfs23;
    RFloat3D *tfsfs11, *tfsfs22, *tfsfs33, *tfsfs12, *tfsfs13, *tfsfs23;
    RFloat3D *rLijMij, *rMijMij;
    RFloat3D *tfsmod;
    RFloat3D *rLijMijavg, *rMijMijavg;
    RFloat3D *cdelta2;

    RFloat1D *fs11, *fs22, *fs33, *fs12, *fs13, *fs23;
    RFloat1D *fsfs11, *fsfs22, *fsfs33, *fsfs12, *fsfs13, *fsfs23;
    RFloat1D *fu1, *fu2, *fu3;
    RFloat1D *fuifuj11, *fuifuj22, *fuifuj33, *fuifuj12, *fuifuj13, *fuifuj23;

    tfu1 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfu2 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfu3 = new RFloat3D(III, JJJ, KKK, fortranArray);

    tfuifuj11 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfuifuj22 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfuifuj33 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfuifuj12 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfuifuj13 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfuifuj23 = new RFloat3D(III, JJJ, KKK, fortranArray);

    tfs11 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfs22 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfs33 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfs12 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfs13 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfs23 = new RFloat3D(III, JJJ, KKK, fortranArray);

    tfsfs11 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfsfs22 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfsfs33 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfsfs12 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfsfs13 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfsfs23 = new RFloat3D(III, JJJ, KKK, fortranArray);

    rLijMij = new RFloat3D(III, JJJ, KKK, fortranArray);
    rMijMij = new RFloat3D(III, JJJ, KKK, fortranArray);

    tfsmod = new RFloat3D(III, JJJ, KKK, fortranArray);

    rLijMijavg = new RFloat3D(III, JJJ, KKK, fortranArray);
    rMijMijavg = new RFloat3D(III, JJJ, KKK, fortranArray);

    cdelta2 = new RFloat3D(III, JJJ, KKK, fortranArray);

    Range IJKMAX (0, MAX(MAX(ni, nj), nk));

    fs11 = new RFloat1D(IJKMAX, fortranArray);
    fs22 = new RFloat1D(IJKMAX, fortranArray);
    fs33 = new RFloat1D(IJKMAX, fortranArray);
    fs12 = new RFloat1D(IJKMAX, fortranArray);
    fs13 = new RFloat1D(IJKMAX, fortranArray);
    fs23 = new RFloat1D(IJKMAX, fortranArray);
    
    fsfs11 = new RFloat1D(IJKMAX, fortranArray);
    fsfs22 = new RFloat1D(IJKMAX, fortranArray);
    fsfs33 = new RFloat1D(IJKMAX, fortranArray);
    fsfs12 = new RFloat1D(IJKMAX, fortranArray);
    fsfs13 = new RFloat1D(IJKMAX, fortranArray);
    fsfs23 = new RFloat1D(IJKMAX, fortranArray);

    fu1 = new RFloat1D(IJKMAX, fortranArray);
    fu2 = new RFloat1D(IJKMAX, fortranArray);
    fu3 = new RFloat1D(IJKMAX, fortranArray);

    fuifuj11 = new RFloat1D(IJKMAX, fortranArray);
    fuifuj22 = new RFloat1D(IJKMAX, fortranArray);
    fuifuj33 = new RFloat1D(IJKMAX, fortranArray);
    fuifuj12 = new RFloat1D(IJKMAX, fortranArray);
    fuifuj13 = new RFloat1D(IJKMAX, fortranArray);
    fuifuj23 = new RFloat1D(IJKMAX, fortranArray);

//
// prepare variables to be filtered
//
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {

                u1 = q(i, j, k, IU);    // q(j,k,i,2:4) = F(u(1:3)) --- IMPLICIT Filter!!!!
                u2 = q(i, j, k, IV);
                u3 = q(i, j, k, IW);

                dudx = gradux(i, j, k);
                dudy = graduy(i, j, k);
                dudz = graduz(i, j, k);

                dvdx = gradvx(i, j, k);
                dvdy = gradvy(i, j, k);
                dvdz = gradvz(i, j, k);

                dwdx = gradwx(i, j, k);
                dwdy = gradwy(i, j, k);
                dwdz = gradwz(i, j, k);

                s11 = dudx;
                s22 = dvdy;
                s33 = dwdz;
                s12 = half * (dudy + dvdx);
                s13 = half * (dudz + dwdx);
                s23 = half * (dvdz + dwdy);

                sij2 = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                smod = sqrt(2.0 * sij2);
                //--- F(ui)
                (*tfu1)(i, j, k) = u1;
                (*tfu2)(i, j, k) = u2;
                (*tfu3)(i, j, k) = u3;
                //--- F(ui)F(uj)
                (*tfuifuj11)(i, j, k) = u1 * u1;
                (*tfuifuj22)(i, j, k) = u2 * u2;
                (*tfuifuj33)(i, j, k) = u3 * u3;
                (*tfuifuj12)(i, j, k) = u1 * u2;
                (*tfuifuj13)(i, j, k) = u1 * u3;
                (*tfuifuj23)(i, j, k) = u2 * u3;
                //--- F(Sij)
                (*tfs11)(i, j, k) = s11;
                (*tfs22)(i, j, k) = s22;
                (*tfs33)(i, j, k) = s33;
                (*tfs12)(i, j, k) = s12;
                (*tfs13)(i, j, k) = s13;
                (*tfs23)(i, j, k) = s23;
                //--- F(|S|)
                (*tfsmod)(i, j, k) = smod;
                //--- F(|S|)F(Sij)
                (*tfsfs11)(i, j, k) = smod * s11;
                (*tfsfs22)(i, j, k) = smod * s22;
                (*tfsfs33)(i, j, k) = smod * s33;
                (*tfsfs12)(i, j, k) = smod * s12;
                (*tfsfs13)(i, j, k) = smod * s13;
                (*tfsfs23)(i, j, k) = smod * s23;
            }
        }
    }

    RDouble a = (*tfs11)(1, 1, 1);
    RDouble b = (*tfs22)(1, 1, 1);
    RDouble c = (*tfs33)(1, 1, 1);
    RDouble d = (*tfs12)(1, 1, 1);
    RDouble e = (*tfs13)(1, 1, 1);
    RDouble f = (*tfs23)(1, 1, 1);
    RDouble g = (*tfsfs11)(1, 1, 1);
    RDouble h = (*tfsfs22)(1, 1, 1);
    RDouble p = (*tfsfs33)(1, 1, 1);
    RDouble r = (*tfsfs12)(1, 1, 1);
    RDouble s = (*tfsfs13)(1, 1, 1);
//
// get filtered variables in I direction
//
    if (filterDirection[0] == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    //--- F(ui)
                    (*fu1)(i) = (*tfu1)(i, j, k);
                    (*fu2)(i) = (*tfu2)(i, j, k);
                    (*fu3)(i) = (*tfu3)(i, j, k);
                    //--- F(ui)F(uj)
                    (*fuifuj11)(i) = (*tfuifuj11)(i, j, k);
                    (*fuifuj22)(i) = (*tfuifuj22)(i, j, k);
                    (*fuifuj33)(i) = (*tfuifuj33)(i, j, k);
                    (*fuifuj12)(i) = (*tfuifuj12)(i, j, k);
                    (*fuifuj13)(i) = (*tfuifuj13)(i, j, k);
                    (*fuifuj23)(i) = (*tfuifuj23)(i, j, k);
                    //--- F(Sij)
                    (*fs11)(i) = (*tfs11)(i, j, k);
                    (*fs22)(i) = (*tfs22)(i, j, k);
                    (*fs33)(i) = (*tfs33)(i, j, k);
                    (*fs12)(i) = (*tfs12)(i, j, k);
                    (*fs13)(i) = (*tfs13)(i, j, k);
                    (*fs23)(i) = (*tfs23)(i, j, k);
                    //--- F(|S|)F(Sij)
                    (*fsfs11)(i) = (*tfsfs11)(i, j, k);
                    (*fsfs22)(i) = (*tfsfs22)(i, j, k);
                    (*fsfs33)(i) = (*tfsfs33)(i, j, k);
                    (*fsfs12)(i) = (*tfsfs12)(i, j, k);
                    (*fsfs13)(i) = (*tfsfs13)(i, j, k);
                    (*fsfs23)(i) = (*tfsfs23)(i, j, k);
                }
                //
                //--- left boundary:
                int i0 = 0;

                u1 = q(i0, j, k, IU);
                u2 = q(i0, j, k, IV);
                u3 = q(i0, j, k, IW);
                //--- F(ui)
                (*fu1)(i0) = u1;
                (*fu2)(i0) = u2;
                (*fu3)(i0) = u3;
                //--- F(ui)F(uj)
                (*fuifuj11)(i0) = u1 * u1;
                (*fuifuj22)(i0) = u2 * u2;
                (*fuifuj33)(i0) = u3 * u3;
                (*fuifuj12)(i0) = u1 * u2;
                (*fuifuj13)(i0) = u1 * u3;
                (*fuifuj23)(i0) = u2 * u3;
                //--- F(Sij) 
                (*fs11)(i0) = (*fs11)(1);
                (*fs22)(i0) = (*fs22)(1);
                (*fs33)(i0) = (*fs33)(1);
                (*fs12)(i0) = (*fs12)(1);
                (*fs13)(i0) = (*fs13)(1);
                (*fs23)(i0) = (*fs23)(1);
                //--- F(|S|)F(Sij)
                (*fsfs11)(i0) = (*fsfs11)(1);
                (*fsfs22)(i0) = (*fsfs22)(1);
                (*fsfs33)(i0) = (*fsfs33)(1);
                (*fsfs12)(i0) = (*fsfs12)(1);
                (*fsfs13)(i0) = (*fsfs13)(1);
                (*fsfs23)(i0) = (*fsfs23)(1);
                //
                //--- right boundary:
                u1 = q(ni, j, k, IU);
                u2 = q(ni, j, k, IV);
                u3 = q(ni, j, k, IW);
                //--- F(ui)
                (*fu1)(ni) = u1;
                (*fu2)(ni) = u2;
                (*fu3)(ni) = u3;
                //--- F(ui)F(uj)
                (*fuifuj11)(ni) = u1 * u1;
                (*fuifuj22)(ni) = u2 * u2;
                (*fuifuj33)(ni) = u3 * u3;
                (*fuifuj12)(ni) = u1 * u2;
                (*fuifuj13)(ni) = u1 * u3;
                (*fuifuj23)(ni) = u2 * u3;
                //--- F(Sij) 
                (*fs11)(ni) = (*fs11)(nis1);
                (*fs22)(ni) = (*fs22)(nis1);
                (*fs33)(ni) = (*fs33)(nis1);
                (*fs12)(ni) = (*fs12)(nis1);
                (*fs13)(ni) = (*fs13)(nis1);
                (*fs23)(ni) = (*fs23)(nis1);
                //---F(|S|)F(Sij)
                (*fsfs11)(ni) = (*fsfs11)(nis1);
                (*fsfs22)(ni) = (*fsfs22)(nis1);
                (*fsfs33)(ni) = (*fsfs33)(nis1);
                (*fsfs12)(ni) = (*fsfs12)(nis1);
                (*fsfs13)(ni) = (*fsfs13)(nis1);
                (*fsfs23)(ni) = (*fsfs23)(nis1);
                //filter
                filter(ni, fu1, fu2, fu3, fuifuj11, fuifuj22, fuifuj33, fuifuj12, fuifuj13, fuifuj23, fs11, fs22, fs33, fs12, fs13, fs23, fsfs11, fsfs22, fsfs33, fsfs12, fsfs13, fsfs23);

                for (int i = ist; i <= ied; ++ i)
                {
                    //--- fill TF(ui)
                    (*tfu1)(i, j, k) = (*fu1)(i);
                    (*tfu2)(i, j, k) = (*fu2)(i);
                    (*tfu3)(i, j, k) = (*fu3)(i);
                    //--- fill T( F(ui)F(uj) )
                    (*tfuifuj11)(i, j, k) = (*fuifuj11)(i);
                    (*tfuifuj22)(i, j, k) = (*fuifuj22)(i);
                    (*tfuifuj33)(i, j, k) = (*fuifuj33)(i);
                    (*tfuifuj12)(i, j, k) = (*fuifuj12)(i);
                    (*tfuifuj13)(i, j, k) = (*fuifuj13)(i);
                    (*tfuifuj23)(i, j, k) = (*fuifuj23)(i);
                    //--- get T( F(Sij) )
                    (*tfs11)(i, j, k) = (*fs11)(i);
                    (*tfs22)(i, j, k) = (*fs22)(i);
                    (*tfs33)(i, j, k) = (*fs33)(i);
                    (*tfs12)(i, j, k) = (*fs12)(i);
                    (*tfs13)(i, j, k) = (*fs13)(i);
                    (*tfs23)(i, j, k) = (*fs23)(i);
                    //--- get T( F(|S|)F(Sij) )
                    (*tfsfs11)(i, j, k) = (*fsfs11)(i);
                    (*tfsfs22)(i, j, k) = (*fsfs22)(i);
                    (*tfsfs33)(i, j, k) = (*fsfs33)(i);
                    (*tfsfs12)(i, j, k) = (*fsfs12)(i);
                    (*tfsfs13)(i, j, k) = (*fsfs13)(i);
                    (*tfsfs23)(i, j, k) = (*fsfs23)(i);
                }
            }
        }
    } //<--- filterDirection[0] == 1

    a = (*tfs11)(1, 1, 1);
    b = (*tfs22)(1, 1, 1);
    c = (*tfs33)(1, 1, 1);
    d = (*tfs12)(1, 1, 1);
    e = (*tfs13)(1, 1, 1);
    f = (*tfs23)(1, 1, 1);
    g = (*tfsfs11)(1, 1, 1);
    h = (*tfsfs22)(1, 1, 1);
    p = (*tfsfs33)(1, 1, 1);
    r = (*tfsfs12)(1, 1, 1);
    s = (*tfsfs13)(1, 1, 1);
//
// get filtered variables in J direction
//
    if (filterDirection[1] == 1) 
    {
        for (int i = ist; i <= ied; ++ i)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    //--- F(ui)
                    (*fu1)(j) = (*tfu1)(i, j, k);
                    (*fu2)(j) = (*tfu2)(i, j, k);
                    (*fu3)(j) = (*tfu3)(i, j, k);
                    //--- F(ui)F(uj)
                    (*fuifuj11)(j) = (*tfuifuj11)(i, j, k);
                    (*fuifuj22)(j) = (*tfuifuj22)(i, j, k);
                    (*fuifuj33)(j) = (*tfuifuj33)(i, j, k);
                    (*fuifuj12)(j) = (*tfuifuj12)(i, j, k);
                    (*fuifuj13)(j) = (*tfuifuj13)(i, j, k);
                    (*fuifuj23)(j) = (*tfuifuj23)(i, j, k);
                    //--- F(Sij)
                    (*fs11)(j) = (*tfs11)(i, j, k);
                    (*fs22)(j) = (*tfs22)(i, j, k);
                    (*fs33)(j) = (*tfs33)(i, j, k);
                    (*fs12)(j) = (*tfs12)(i, j, k);
                    (*fs13)(j) = (*tfs13)(i, j, k);
                    (*fs23)(j) = (*tfs23)(i, j, k);
                    //--- F(|S|)F(Sij)
                    (*fsfs11)(j) = (*tfsfs11)(i, j, k);
                    (*fsfs22)(j) = (*tfsfs22)(i, j, k);
                    (*fsfs33)(j) = (*tfsfs33)(i, j, k);
                    (*fsfs12)(j) = (*tfsfs12)(i, j, k);
                    (*fsfs13)(j) = (*tfsfs13)(i, j, k);
                    (*fsfs23)(j) = (*tfsfs23)(i, j, k);
                }
                //
                //--- left boundary:
                int j0 = 0;

                u1 = q(i, j0, k, IU);
                u2 = q(i, j0, k, IV);
                u3 = q(i, j0, k, IW);
                //--- F(ui)
                (*fu1)(j0) = u1;
                (*fu2)(j0) = u2;
                (*fu3)(j0) = u3;
                //--- F(ui)F(uj)
                (*fuifuj11)(j0) = u1 * u1;
                (*fuifuj22)(j0) = u2 * u2;
                (*fuifuj33)(j0) = u3 * u3;
                (*fuifuj12)(j0) = u1 * u2;
                (*fuifuj13)(j0) = u1 * u3;
                (*fuifuj23)(j0) = u2 * u3;
                //--- F(Sij) 
                (*fs11)(j0) = (*fs11)(1);
                (*fs22)(j0) = (*fs22)(1);
                (*fs33)(j0) = (*fs33)(1);
                (*fs12)(j0) = (*fs12)(1);
                (*fs13)(j0) = (*fs13)(1);
                (*fs23)(j0) = (*fs23)(1);
                //---F(|S|)F(Sij)
                (*fsfs11)(j0) = (*fsfs11)(1);
                (*fsfs22)(j0) = (*fsfs22)(1);
                (*fsfs33)(j0) = (*fsfs33)(1);
                (*fsfs12)(j0) = (*fsfs12)(1);
                (*fsfs13)(j0) = (*fsfs13)(1);
                (*fsfs23)(j0) = (*fsfs23)(1);
                //
                //--- right boundary:
                u1 = q(i, nj, k, IU);
                u2 = q(i, nj, k, IV);
                u3 = q(i, nj, k, IW);
                //--- F(ui)
                (*fu1)(nj) = u1;
                (*fu2)(nj) = u2;
                (*fu3)(nj) = u3;
                //--- F(ui)F(uj)
                (*fuifuj11)(nj) = u1 * u1;
                (*fuifuj22)(nj) = u2 * u2;
                (*fuifuj33)(nj) = u3 * u3;
                (*fuifuj12)(nj) = u1 * u2;
                (*fuifuj13)(nj) = u1 * u3;
                (*fuifuj23)(nj) = u2 * u3;
                //--- F(Sij)
                (*fs11)(nj) = (*fs11)(njs1);
                (*fs22)(nj) = (*fs22)(njs1);
                (*fs33)(nj) = (*fs33)(njs1);
                (*fs12)(nj) = (*fs12)(njs1);
                (*fs13)(nj) = (*fs13)(njs1);
                (*fs23)(nj) = (*fs23)(njs1);
                //---F(|S|)F(Sij)
                (*fsfs11)(nj) = (*fsfs11)(njs1);
                (*fsfs22)(nj) = (*fsfs22)(njs1);
                (*fsfs33)(nj) = (*fsfs33)(njs1);
                (*fsfs12)(nj) = (*fsfs12)(njs1);
                (*fsfs13)(nj) = (*fsfs13)(njs1);
                (*fsfs23)(nj) = (*fsfs23)(njs1);        
                //filter
                filter(nj, fu1, fu2, fu3, fuifuj11, fuifuj22, fuifuj33, fuifuj12, fuifuj13, fuifuj23, fs11, fs22, fs33, fs12, fs13, fs23, fsfs11, fsfs22, fsfs33, fsfs12, fsfs13, fsfs23);

                for (int j = jst; j <= jed; ++ j)
                {
                    //--- fill TF(ui)
                    (*tfu1)(i,j,k) = (*fu1)(j);
                    (*tfu2)(i,j,k) = (*fu2)(j);
                    (*tfu3)(i,j,k) = (*fu3)(j);
                    //--- fill T( F(ui)F(uj) )
                    (*tfuifuj11)(i,j,k) = (*fuifuj11)(j);
                    (*tfuifuj22)(i,j,k) = (*fuifuj22)(j);
                    (*tfuifuj33)(i,j,k) = (*fuifuj33)(j);
                    (*tfuifuj12)(i,j,k) = (*fuifuj12)(j);
                    (*tfuifuj13)(i,j,k) = (*fuifuj13)(j);
                    (*tfuifuj23)(i,j,k) = (*fuifuj23)(j);
                    //--- get T( F(Sij) )
                    (*tfs11)(i,j,k) = (*fs11)(j);
                    (*tfs22)(i,j,k) = (*fs22)(j);
                    (*tfs33)(i,j,k) = (*fs33)(j);
                    (*tfs12)(i,j,k) = (*fs12)(j);
                    (*tfs13)(i,j,k) = (*fs13)(j);
                    (*tfs23)(i,j,k) = (*fs23)(j);
                    //--- get T( F(|S|)F(Sij) )
                    (*tfsfs11)(i,j,k) = (*fsfs11)(j);
                    (*tfsfs22)(i,j,k) = (*fsfs22)(j);
                    (*tfsfs33)(i,j,k) = (*fsfs33)(j);
                    (*tfsfs12)(i,j,k) = (*fsfs12)(j);
                    (*tfsfs13)(i,j,k) = (*fsfs13)(j);
                    (*tfsfs23)(i,j,k) = (*fsfs23)(j);
                }
            }
        }

    } //<--- filterDirection[1] == 0

    a = (*tfs11)(1, 1, 1);
    b = (*tfs22)(1, 1, 1);
    c = (*tfs33)(1, 1, 1);
    d = (*tfs12)(1, 1, 1);
    e = (*tfs13)(1, 1, 1);
    f = (*tfs23)(1, 1, 1);
    g = (*tfsfs11)(1, 1, 1);
    h = (*tfsfs22)(1, 1, 1);
    p = (*tfsfs33)(1, 1, 1);
    r = (*tfsfs12)(1, 1, 1);
    s = (*tfsfs13)(1, 1, 1);
//
// get filtered variables in K direction
//
    if (filterDirection[2] == 1) 
    {
        for (int i = ist; i <= ied; ++ i)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    //--- F(ui)
                    (*fu1)(k) = (*tfu1)(i,j,k);
                    (*fu2)(k) = (*tfu2)(i,j,k);
                    (*fu3)(k) = (*tfu3)(i,j,k);
                    //--- F(ui)F(uj)
                    (*fuifuj11)(k) = (*tfuifuj11)(i,j,k);
                    (*fuifuj22)(k) = (*tfuifuj22)(i,j,k);
                    (*fuifuj33)(k) = (*tfuifuj33)(i,j,k);
                    (*fuifuj12)(k) = (*tfuifuj12)(i,j,k);
                    (*fuifuj13)(k) = (*tfuifuj13)(i,j,k);
                    (*fuifuj23)(k) = (*tfuifuj23)(i,j,k);
                    //--- F(Sij)
                    (*fs11)(k) = (*tfs11)(i,j,k);
                    (*fs22)(k) = (*tfs22)(i,j,k);
                    (*fs33)(k) = (*tfs33)(i,j,k);
                    (*fs12)(k) = (*tfs12)(i,j,k);
                    (*fs13)(k) = (*tfs13)(i,j,k);
                    (*fs23)(k) = (*tfs23)(i,j,k);
                    //--- F(|S|)F(Sij)
                    (*fsfs11)(k) = (*tfsfs11)(i,j,k);
                    (*fsfs22)(k) = (*tfsfs22)(i,j,k);
                    (*fsfs33)(k) = (*tfsfs33)(i,j,k);
                    (*fsfs12)(k) = (*tfsfs12)(i,j,k);
                    (*fsfs13)(k) = (*tfsfs13)(i,j,k);
                    (*fsfs23)(k) = (*tfsfs23)(i,j,k);
                }
                //
                //--- left boundary:
                int k0 = 0;

                u1 = q(i,j,k0,IU);
                u2 = q(i,j,k0,IV);
                u3 = q(i,j,k0,IW);
                //--- F(ui)
                (*fu1)(k0) = u1;
                (*fu2)(k0) = u2;
                (*fu3)(k0) = u3;
                //--- F(ui)F(uj)
                (*fuifuj11)(k0) = u1*u1;
                (*fuifuj22)(k0) = u2*u2;
                (*fuifuj33)(k0) = u3*u3;
                (*fuifuj12)(k0) = u1*u2;
                (*fuifuj13)(k0) = u1*u3;
                (*fuifuj23)(k0) = u2*u3;
                //--- F(Sij) 
                (*fs11)(k0) = (*fs11)(1);
                (*fs22)(k0) = (*fs22)(1);
                (*fs33)(k0) = (*fs33)(1);
                (*fs12)(k0) = (*fs12)(1);
                (*fs13)(k0) = (*fs13)(1);
                (*fs23)(k0) = (*fs23)(1);
                //--- F(|S|)F(Sij)
                (*fsfs11)(k0) = (*fsfs11)(1);
                (*fsfs22)(k0) = (*fsfs22)(1);
                (*fsfs33)(k0) = (*fsfs33)(1);
                (*fsfs12)(k0) = (*fsfs12)(1);
                (*fsfs13)(k0) = (*fsfs13)(1);
                (*fsfs23)(k0) = (*fsfs23)(1);
                //
                //--- right boundary:
                u1 = q(i,j,nk,IU);
                u2 = q(i,j,nk,IV);
                u3 = q(i,j,nk,IW);
                //--- F(ui)
                (*fu1)(nk) = u1;
                (*fu2)(nk) = u2;
                (*fu3)(nk) = u3;
                //--- F(ui)F(uj)
                (*fuifuj11)(nk) = u1*u1;
                (*fuifuj22)(nk) = u2*u2;
                (*fuifuj33)(nk) = u3*u3;
                (*fuifuj12)(nk) = u1*u2;
                (*fuifuj13)(nk) = u1*u3;
                (*fuifuj23)(nk) = u2*u3;
                //--- F(Sij)
                (*fs11)(nk) = (*fs11)(nks1);
                (*fs22)(nk) = (*fs22)(nks1);
                (*fs33)(nk) = (*fs33)(nks1);
                (*fs12)(nk) = (*fs12)(nks1);
                (*fs13)(nk) = (*fs13)(nks1);
                (*fs23)(nk) = (*fs23)(nks1);
                //--- F(|S|)F(Sij)
                (*fsfs11)(nk) = (*fsfs11)(nks1);
                (*fsfs22)(nk) = (*fsfs22)(nks1);
                (*fsfs33)(nk) = (*fsfs33)(nks1);
                (*fsfs12)(nk) = (*fsfs12)(nks1);
                (*fsfs13)(nk) = (*fsfs13)(nks1);
                (*fsfs23)(nk) = (*fsfs23)(nks1);        
                //filter
                filter(nk,fu1,fu2,fu3,fuifuj11,fuifuj22,fuifuj33,fuifuj12,fuifuj13,fuifuj23,fs11,fs22,fs33,fs12,fs13,fs23,fsfs11,fsfs22,fsfs33,fsfs12,fsfs13,fsfs23);
        
                for (int k = kst; k <= ked; ++ k)
                {
                    //--- fill TF(ui)
                    (*tfu1)(i,j,k) = (*fu1)(k);
                    (*tfu2)(i,j,k) = (*fu2)(k);
                    (*tfu3)(i,j,k) = (*fu3)(k);
                    //--- fill T( F(ui)F(uj) )
                    (*tfuifuj11)(i,j,k) = (*fuifuj11)(k);
                    (*tfuifuj22)(i,j,k) = (*fuifuj22)(k);
                    (*tfuifuj33)(i,j,k) = (*fuifuj33)(k);
                    (*tfuifuj12)(i,j,k) = (*fuifuj12)(k);
                    (*tfuifuj13)(i,j,k) = (*fuifuj13)(k);
                    (*tfuifuj23)(i,j,k) = (*fuifuj23)(k);
                    //--- get T( F(Sij) )
                    (*tfs11)(i,j,k) = (*fs11)(k);
                    (*tfs22)(i,j,k) = (*fs22)(k);
                    (*tfs33)(i,j,k) = (*fs33)(k);
                    (*tfs12)(i,j,k) = (*fs12)(k);
                    (*tfs13)(i,j,k) = (*fs13)(k);
                    (*tfs23)(i,j,k) = (*fs23)(k);
                    //--- get T( F(|S|)F(Sij) )
                    (*tfsfs11)(i,j,k) = (*fsfs11)(k);
                    (*tfsfs22)(i,j,k) = (*fsfs22)(k);
                    (*tfsfs33)(i,j,k) = (*fsfs33)(k);
                    (*tfsfs12)(i,j,k) = (*fsfs12)(k);
                    (*tfsfs13)(i,j,k) = (*fsfs13)(k);
                    (*tfsfs23)(i,j,k) = (*fsfs23)(k);
                }
            }
        }

    } //<--- filterDirection[2] == 1
//
// get LijMij, MijMij
    for (int i = ist; i <= ied; ++ i) 
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                //--- Lij = T( F(ui)F(uj) ) - TF(ui)*TF(uj)
                rL11 = (*tfuifuj11)(i,j,k) - (*tfu1)(i,j,k) * (*tfu1)(i,j,k);
                rL22 = (*tfuifuj22)(i,j,k) - (*tfu2)(i,j,k) * (*tfu2)(i,j,k);
                rL33 = (*tfuifuj33)(i,j,k) - (*tfu3)(i,j,k) * (*tfu3)(i,j,k);
                rL12 = (*tfuifuj12)(i,j,k) - (*tfu1)(i,j,k) * (*tfu2)(i,j,k);
                rL13 = (*tfuifuj13)(i,j,k) - (*tfu1)(i,j,k) * (*tfu3)(i,j,k);
                rL23 = (*tfuifuj23)(i,j,k) - (*tfu2)(i,j,k) * (*tfu3)(i,j,k);
                //--- Mij = ( TF(Detla)/F(Delta) )^2*|TF(S)|*T( F(Sij) ) -  T( F(|S|)F(Sij) )
                //--- TF(|S|)
                tfsijtfsij = (*tfs11)(i,j,k) * (*tfs11)(i,j,k) 
                           + (*tfs22)(i,j,k) * (*tfs22)(i,j,k)
                           + (*tfs33)(i,j,k) * (*tfs33)(i,j,k)
                           + 2.0*(  (*tfs12)(i,j,k) * (*tfs12)(i,j,k)
                                  + (*tfs13)(i,j,k) * (*tfs13)(i,j,k)
                                  + (*tfs23)(i,j,k) * (*tfs23)(i,j,k) );
                tfss = sqrt(2.0*tfsijtfsij);

                rM11 = rfwr * rfwr * tfss * (*tfs11)(i,j,k) - (*tfsfs11)(i,j,k);
                rM22 = rfwr * rfwr * tfss * (*tfs22)(i,j,k) - (*tfsfs22)(i,j,k);
                rM33 = rfwr * rfwr * tfss * (*tfs33)(i,j,k) - (*tfsfs33)(i,j,k);
                rM12 = rfwr * rfwr * tfss * (*tfs12)(i,j,k) - (*tfsfs12)(i,j,k);
                rM13 = rfwr * rfwr * tfss * (*tfs13)(i,j,k) - (*tfsfs13)(i,j,k);
                rM23 = rfwr * rfwr * tfss * (*tfs23)(i,j,k) - (*tfsfs23)(i,j,k);

                (*rLijMij)(i,j,k) = rL11*rM11 +rL22*rM22 +rL33*rM33
                                  + 2.0*( rL12*rM12 +rL13*rM13 +rL23*rM23 );

                (*rMijMij)(i,j,k) = rM11*rM11 +rM22*rM22 +rM33*rM33
                                  + 2.0*( rM12*rM12 +rM13*rM13 +rM23*rM23 );
            }
        }
    }

    a = (*tfs11)(1, 1, 1);
    b = (*tfsfs11)(1, 1, 1);
//
// get <LijMij>, <MijMij> in I direction
    
    if ((averageDirection[0] == 0) && (averageDirection[1] == 0) && (averageDirection[2] == 0))
    {
        int averageWidth = parameters->GetAverageWidth();

        int naw = averageWidth;
        //--- Boundary points and corner points were averaged only with points inside this block!
        //    This should be OK, since the foregoing filter operations were executed only on points in this block also!
        for (int i = ist; i <= ied; ++ i)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    RDouble rLMavg = 0.0;
                    RDouble rMMavg = 0.0;
                    int nnavg = 0;
                    for (int idi = -naw; idi <= naw; ++ idi)
                    {
                        for (int idj = -naw; j <= naw; ++ idj)
                        {
                            for (int idk = -naw; k <= naw; ++ idk)
                            {
                                int ii = i + idi;
                                int jj = j + idj;
                                int kk = k + idk;
                                bool Lvalid = (ii>=1 && ii<=nis1) && (jj>=1 && jj<=njs1) && (kk>=1 && kk<=nks1);
                                if (Lvalid) 
                                {
                                    nnavg = nnavg + 1;
                                    rLMavg = rLMavg + (*rLijMij)(jj,kk,ii);
                                    rMMavg = rMMavg + (*rMijMij)(jj,kk,ii);
                                }
                            }
                        }
                    }
                    (*rLijMijavg)(i,j,k) = rLMavg/nnavg;
                    (*rMijMijavg)(i,j,k) = rMMavg/nnavg;
                }
            }
        }

        for (int i = ist; i <= ied; ++ i)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    (*rLijMij)(i,j,k) = (*rLijMijavg)(i,j,k);
                    (*rMijMij)(i,j,k) = (*rMijMijavg)(i,j,k);
                }
            }
        }
    }
    else
    {
        if (averageDirection[0] != 0)
        // get <LijMij>, <MijMij> in I direction
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    RDouble rLMavg = 0.0;
                    RDouble rMMavg = 0.0;
                    for (int i = ist; i <= ied; ++ i)
                    {
                        rLMavg = rLMavg + (*rLijMij)(i,j,k);
                        rMMavg = rMMavg + (*rMijMij)(i,j,k);
                    }
                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*rLijMij)(i,j,k) = rLMavg/nis1;
                        (*rMijMij)(i,j,k) = rMMavg/nis1;
                    }
                }
            }
        }

        RDouble lm = (*rLijMij)(1, 1, 1);
        RDouble mm = (*rMijMij)(1, 1, 1);

        if (averageDirection[1] != 0)
        // get <LijMij>, <MijMij> in J direction
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble rLMavg = 0.0;
                    RDouble rMMavg = 0.0;
                    for (int j = ist; j <= jed; ++ j)
                    {
                        RDouble a = (*rLijMij)(i,j,k);
                        RDouble b = (*rMijMij)(i,j,k);
                        rLMavg = rLMavg + (*rLijMij)(i,j,k);
                        rMMavg = rMMavg + (*rMijMij)(i,j,k);
                    }
                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*rLijMij)(i,j,k) = rLMavg/njs1;
                        (*rMijMij)(i,j,k) = rMMavg/njs1;
                    }
                }
            }
        }

        lm = (*rLijMij)(1, 1, 1);
        mm = (*rMijMij)(1, 1, 1);

        if (averageDirection[2] != 0)
        // get <LijMij>, <MijMij> in K direction
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    RDouble rLMavg = 0.0;
                    RDouble rMMavg = 0.0;
                    for (int k = kst; k <= ked; ++ k)
                    {
                        rLMavg = rLMavg + (*rLijMij)(i,j,k);
                        rMMavg = rMMavg + (*rMijMij)(i,j,k);
                    }
                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*rLijMij)(i,j,k) = rLMavg/nks1;
                        (*rMijMij)(i,j,k) = rMMavg/nks1;
                    }
                }
            }
        }
    }

//
// get cdelta2, get vist(i,j,k)

    int nNegativeCr = 0;
    int nNegativeMu = 0;
    int nLowerBoundCr = 0;
    int nUpperBoundCr = 0;

    for (int i = ist; i <= ied; ++ i)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                RDouble lijmij = (*rLijMij)(i,j,k);
                RDouble mijmij = (*rMijMij)(i,j,k);
                RDouble cr = -0.5 * lijmij / mijmij;

                (*cdelta2)(i,j,k) = cr;

                u1 = q(i,j,k,IU);    // q(j,k,i,2:4) = F( u(1:3) ) --- IMPLICIT Filter!!!!
                u2 = q(i,j,k,IV);
                u3 = q(i,j,k,IW);

                dudx = gradux(i,j,k);
                dudy = graduy(i,j,k);
                dudz = graduz(i,j,k);

                dvdx = gradvx(i,j,k);
                dvdy = gradvy(i,j,k);
                dvdz = gradvz(i,j,k);

                dwdx = gradwx(i,j,k);
                dwdy = gradwy(i,j,k);
                dwdz = gradwz(i,j,k);

                s11 = dudx;
                s22 = dvdy;
                s33 = dwdz;
                s12 = half * (dudy + dvdx);
                s13 = half * (dudz + dwdx);
                s23 = half * (dvdz + dwdy);

                sij2 = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                smod = sqrt(2.0*sij2);
                rho = q(i,j,k,IR);
                mul = visl(i,j,k);

                mut = rho * refReNumber * (*cdelta2)(i,j,k) * smod;
                //--- smagorinsky nut

                if (cr < 0.0)
                {
                    nNegativeCr = nNegativeCr + 1;
                }

                if (cr < -0.01)
                {
                    //cr = -0.01;
                    nLowerBoundCr = nLowerBoundCr + 1;
                }

                if (cr > 0.5)
                {
                    //cr = 0.5;
                    nUpperBoundCr = nUpperBoundCr + 1;
                }

                if (mut < -mul)
                {
                    nNegativeMu = nNegativeMu + 1;
                }

                if (turbViscousCutType == 1)
                {
                    mut = MAX(mut - mul, 0.0);
                }
                else if (turbViscousCutType == 2)
                {
                    mut = MAX(mut, 0.0);
                }

                vist(i,j,k) = mut;
            }
        }
    }

    grid->UpdateData("nNegativeCr", &nNegativeCr, PHINT, 1);
    grid->UpdateData("nNegativeMu", &nNegativeMu, PHINT, 1);

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");
    int monitorNegativeConstant = GlobalDataBase::GetIntParaFromDB("monitorNegativeConstant");

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorNegativeConstant)
        {
            RDouble percentNegativeCr = nNegativeCr / static_cast<RDouble> ((ni-1) * (nj-1) * (nk-1));
            RDouble percentNegativeMu = nNegativeMu / static_cast<RDouble> ((ni-1) * (nj-1) * (nk-1));
            RDouble percentLowerBoundCr = nLowerBoundCr / static_cast<RDouble> ((ni-1) * (nj-1) * (nk-1));
            RDouble percentUpperBoundCr = nUpperBoundCr / static_cast<RDouble> ((ni-1) * (nj-1) * (nk-1));

            cout << "percentNegativeCr = " << percentNegativeCr << "\n";
            //cout << "nNegativeCi = " << nNegativeCi << "\n";
            //cout << "nNegativePrt = " << nNegativePrt << "\n";
            cout << "percentNegativeMu = " << percentNegativeMu << "\n";
            cout << "percentLowerBoundCr = " << percentLowerBoundCr << "\n";
            cout << "percentUpperBoundCr = " << percentUpperBoundCr << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif

    delete tfu1      ; delete tfu2      ; delete tfu3;
    delete tfuifuj11 ; delete tfuifuj22 ; delete tfuifuj33; delete tfuifuj12; delete tfuifuj13; delete tfuifuj23;
    delete tfs11     ; delete tfs22     ; delete tfs33    ; delete tfs12    ; delete tfs13    ; delete tfs23    ;
    delete tfsfs11   ; delete tfsfs22   ; delete tfsfs33  ; delete tfsfs12  ; delete tfsfs13  ; delete tfsfs23  ;
    delete rLijMij   ; delete rMijMij   ;
    delete tfsmod    ;
    delete rLijMijavg; delete rMijMijavg;
    delete cdelta2   ;

    delete fs11      ; delete fs22      ; delete fs33     ; delete fs12     ; delete fs13     ; delete fs23     ;
    delete fsfs11    ; delete fsfs22    ; delete fsfs33   ; delete fsfs12   ; delete fsfs13   ; delete fsfs23   ;
    delete fu1       ; delete fu2       ; delete fu3      ;
    delete fuifuj11  ; delete fuifuj22  ; delete fuifuj33 ; delete fuifuj12 ; delete fuifuj13 ; delete fuifuj23 ;
}


void LESSolverStruct::DynamicSmagViscosityCompressibleOld(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble refReNumber = parameters->GetRefReNumber();

    int turbViscousCutType = parameters->GetTurbViscousCutType();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    int monitorNegativeConstant = GlobalDataBase::GetIntParaFromDB("monitorNegativeConstant");

    int prtCutType = GlobalDataBase::GetIntParaFromDB("prtCutType");

    RDouble oprt = parameters->GetoPrandtlTurbulence();

    RDouble refGama = parameters->GetRefGama();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);

    RDouble4D &wallFunction = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble5D &gradPrimtiveVarFaceXofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceXofLES"));
    RDouble5D &gradPrimtiveVarFaceYofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceYofLES"));
    RDouble5D &gradPrimtiveVarFaceZofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceZofLES"));

    RDouble5D &gradtx = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceXofLES"));
    RDouble5D &gradty = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceYofLES"));
    RDouble5D &gradtz = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceZofLES"));

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));

    RDouble5D &subgridStress = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridStress"));
    //RDouble4D &subgridEnergy = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("kSGS"));
    RDouble5D &subgridHeatFlux = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridHeatFlux"));

    ComputeAnisotropicConstant(grid);

    ComputeIsotropicConstant(grid);
    //ComputeSubgridScaleEnergy(grid);

    ComputePrandtlNumber(grid);

    double vistmax  = 0.0;
    double vistmin  = 1.0e30;
    double logsmall = log(SMALL);

    int nNegativeCr = 0;
    int nLowerBoundCr = 0;
    int nUpperBoundCr = 0;
    int nNegativeCi = 0;
    int nNegativePrt = 0;
    int nNegativeMu = 0;

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    int imin = 1;
    int jmin = 1;
    int kmin = 1;

    int nDim = GetDim();

    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int ist, ied, jst, jed, kst, ked;
        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble dudx = gradPrimtiveVarFaceXofLES(il, jl, kl, IU, iSurface);
                    RDouble dudy = gradPrimtiveVarFaceYofLES(il, jl, kl, IU, iSurface);
                    RDouble dudz = gradPrimtiveVarFaceZofLES(il, jl, kl, IU, iSurface);

                    RDouble dvdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IV, iSurface);

                    RDouble dwdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IW, iSurface);

                    RDouble dtdx = gradtx(il, jl, kl, ITT, iSurface);
                    RDouble dtdy = gradty(il, jl, kl, ITT, iSurface);
                    RDouble dtdz = gradtz(il, jl, kl, ITT, iSurface);

                    RDouble s11 = dudx;
                    RDouble s22 = dvdy;
                    RDouble s33 = dwdz;
                    RDouble s12 = half * (dudy + dvdx);
                    RDouble s13 = half * (dudz + dwdx);
                    RDouble s23 = half * (dvdz + dwdy);

                    RDouble sijsij = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                    RDouble ss = sqrt(two * sijsij);
                    RDouble skk = s11 + s22 + s33;

                    RDouble rho = 0.5 * (q(i, j, k, IR) + q(il, jl, kl, IR));
                    RDouble mul = 0.5 * (visl(i, j, k) + visl(il, jl, kl));
                    RDouble dampingFunction = wallFunction(il, jl, kl, iSurface);

                    RDouble cr = 0.5 * (anisotropicConstant(il, jl, kl) + anisotropicConstant(i, j, k));

                    if (cr < -0.01)
                    {
                        //cr = -0.01;
                        nLowerBoundCr = nLowerBoundCr + 1;
                    }

                    /*
                    if (cr > 0.5)
                    {
                        cr = 0.5;
                        nUpperBoundCr = nUpperBoundCr + 1;
                    }
                    */

                    RDouble mut = cr * dampingFunction * dampingFunction * rho * ss * refReNumber;

                    RDouble ci = 0.5 * (isotropicConstant(i, j, k) + isotropicConstant(il, jl, kl));
                    //ci = MIN(MAX(ci, 0.0), 0.01);

                    RDouble kSGS = 2.0 * ci * dampingFunction * dampingFunction * rho * ss * ss * refReNumber;

                    oprt = 2.0 / (turbulentPrandtlNumber(il, jl, kl) + turbulentPrandtlNumber(i, j, k));
                    //kSGS = 0.0;

                    RDouble hfxSGS = -mut * specificHeatAtConstantPressure * oprt * dtdx;
                    RDouble hfySGS = -mut * specificHeatAtConstantPressure * oprt * dtdy;
                    RDouble hfzSGS = -mut * specificHeatAtConstantPressure * oprt * dtdz;

                    subgridHeatFlux(il, jl, kl, 1, iSurface) = hfxSGS;
                    subgridHeatFlux(il, jl, kl, 2, iSurface) = hfySGS;
                    subgridHeatFlux(il, jl, kl, 3, iSurface) = hfzSGS;

                    if (mut < -mul)
                    {
                        nNegativeMu = nNegativeMu + 1;
                    }

                    if (turbViscousCutType == 1)
                    {
                        mut = MAX(mut - mul, 0.0);
                    }
                    else if (turbViscousCutType == 2)
                    {
                        mut = MAX(mut, 0.0);
                    }
                    else if (turbViscousCutType == 3)
                    {
                        mut = MAX(mut, -mul);
                    }

                    if (prtCutType == 1)
                    {
                        oprt = MAX(oprt, 0.0);
                    }

                    mut = MIN(eddyViscosityLimit, mut);

                    RDouble txxSGS = -2.0 * mut * (s11 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    RDouble tyySGS = -2.0 * mut * (s22 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    RDouble tzzSGS = -2.0 * mut * (s33 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    RDouble txySGS = -2.0 * mut * s12;
                    RDouble txzSGS = -2.0 * mut * s13;
                    RDouble tyzSGS = -2.0 * mut * s23;

                    subgridStress(il, jl, kl, 1, iSurface) = txxSGS;
                    subgridStress(il, jl, kl, 2, iSurface) = tyySGS;
                    subgridStress(il, jl, kl, 3, iSurface) = tzzSGS;
                    subgridStress(il, jl, kl, 4, iSurface) = txySGS;
                    subgridStress(il, jl, kl, 5, iSurface) = txzSGS;
                    subgridStress(il, jl, kl, 6, iSurface) = tyzSGS;     

                    if (vistmax < mut)
                    {
                        vistmax = mut;
                        imax = i;
                        jmax = j;
                        kmax = k;
                    }

                    if (vistmin > mut)
                    {
                        vistmin = mut;
                        imin = i;
                        jmin = j;
                        kmin = k;
                    }

                    if (cr < 0.0)
                    {
                        nNegativeCr = nNegativeCr + 1;
                    }

                    if (ci < 0.0)
                    {
                        nNegativeCi = nNegativeCi + 1;
                    }

                    if (oprt < 0.0)
                    {
                        nNegativePrt = nNegativePrt + 1;
                    }
                }
            }
        }

    }

    double vist_max = vistmax;
    double vist_min = vistmin;

    grid->UpdateData("vist_max", &vist_max, PHDOUBLE, 1);
    grid->UpdateData("vist_min", &vist_min, PHDOUBLE, 1);
    grid->UpdateData("nNegativeCr", &nNegativeCr, PHINT, 1);
    grid->UpdateData("nNegativeCi", &nNegativeCi, PHINT, 1);
    grid->UpdateData("nNegativePrt", &nNegativePrt, PHINT, 1);
    grid->UpdateData("nNegativeMu", &nNegativeMu, PHINT, 1);
    grid->UpdateData("nLowerBoundCr", &nLowerBoundCr, PHINT, 1);
    grid->UpdateData("nUpperBoundCr", &nUpperBoundCr, PHINT, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorVistmax)
        {
            cout << "vistmax = " << vistmax << " " << imax << " " << jmax << " " << kmax << "\n";
            cout << "vistmin = " << vistmin << " " << imin << " " << jmin << " " << kmin << "\n";
        }

        if (monitorNegativeConstant)
        {
            RDouble percentNegativeCr = nNegativeCr / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentNegativeMu = nNegativeMu / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentLowerBoundCr = nLowerBoundCr / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentUpperBoundCr = nUpperBoundCr / static_cast<RDouble> (ni * nj * nk * nDim);

            cout << "percentNegativeCr = " << percentNegativeCr << "\n";
            //cout << "nNegativeCi = " << nNegativeCi << "\n";
            //cout << "nNegativePrt = " << nNegativePrt << "\n";
            cout << "percentNegativeMu = " << percentNegativeMu << "\n";
            cout << "percentLowerBoundCr = " << percentLowerBoundCr << "\n";
            cout << "percentUpperBoundCr = " << percentUpperBoundCr << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void LESSolverStruct::DynamicSmagViscosityCompressible(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_LESSolverStruct *parameters = GetControlParameters();

    //RDouble refReNumber = parameters->GetRefReNumber();

    int monitorNegativeConstant = GlobalDataBase::GetIntParaFromDB("monitorNegativeConstant");

    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));

    ComputeAnisotropicConstant(grid);

    ComputeIsotropicConstant(grid);
    //ComputeSubgridScaleEnergy(grid);

    ComputePrandtlNumber(grid);

    double vistmax  = 0.0;
    double vistmin  = 1.0e30;
    double logsmall = log(SMALL);

    int nNegativeCr = 0;
    int nLowerBoundCr = 0;
    int nUpperBoundCr = 0;
    int nNegativeCi = 0;
    int nNegativePrt = 0;

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    int imin = 1;
    int jmin = 1;
    int kmin = 1;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble dampingFunction = wallFunction(i, j, k);

                RDouble cr = anisotropicConstant(i, j, k);

                RDouble ci = isotropicConstant(i, j, k);

                RDouble prt = turbulentPrandtlNumber(i, j, k);

                if (cr < -0.01)
                {
                    //cr = -0.01;
                    nLowerBoundCr = nLowerBoundCr + 1;
                }

                if (cr > 0.5)
                {
                    //cr = 0.5;
                    nUpperBoundCr = nUpperBoundCr + 1;
                }

                if (cr < 0.0)
                {
                    nNegativeCr = nNegativeCr + 1;
                }

                if (ci < 0.0)
                {
                    nNegativeCi = nNegativeCi + 1;
                }

                if (prt < 0.0)
                {
                    nNegativePrt = nNegativePrt + 1;
                }

                anisotropicConstant(i, j, k) = cr * dampingFunction * dampingFunction;

                isotropicConstant(i, j, k) = ci * dampingFunction * dampingFunction;

                //turbulentPrandtlNumber(i, j, k) = prt;
            }
        }
    }

    grid->UpdateData("nNegativeCr", &nNegativeCr, PHINT, 1);
    grid->UpdateData("nNegativeCi", &nNegativeCi, PHINT, 1);
    grid->UpdateData("nNegativePrt", &nNegativePrt, PHINT, 1);
    grid->UpdateData("nLowerBoundCr", &nLowerBoundCr, PHINT, 1);
    grid->UpdateData("nUpperBoundCr", &nUpperBoundCr, PHINT, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorNegativeConstant)
        {
            RDouble percentNegativeCr = nNegativeCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentNegativeCi = nNegativeCi / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentNegativePrt = nNegativePrt / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentLowerBoundCr = nLowerBoundCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentUpperBoundCr = nUpperBoundCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));

            cout << "percentNegativeCr = " << percentNegativeCr << "\n";
            cout << "percentNegativeCi = " << nNegativeCi << "\n";
            cout << "percentNegativePrt = " << nNegativePrt << "\n";
            cout << "percentLowerBoundCr = " << percentLowerBoundCr << "\n";
            cout << "percentUpperBoundCr = " << percentUpperBoundCr << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void LESSolverStruct::DynamicSmagViscosityIncompressibleOld(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble refReNumber = parameters->GetRefReNumber();

    int turbViscousCutType = parameters->GetTurbViscousCutType();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    int monitorNegativeConstant = GlobalDataBase::GetIntParaFromDB("monitorNegativeConstant");

    int prtCutType = GlobalDataBase::GetIntParaFromDB("prtCutType");

    RDouble oprt = parameters->GetoPrandtlTurbulence();

    RDouble refGama = parameters->GetRefGama();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);

    RDouble4D &wallFunction = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble5D &gradPrimtiveVarFaceXofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceXofLES"));
    RDouble5D &gradPrimtiveVarFaceYofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceYofLES"));
    RDouble5D &gradPrimtiveVarFaceZofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceZofLES"));

    RDouble5D &gradtx = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceXofLES"));
    RDouble5D &gradty = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceYofLES"));
    RDouble5D &gradtz = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceZofLES"));

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));

    RDouble5D &subgridStress = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridStress"));
    RDouble5D &subgridHeatFlux = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridHeatFlux"));

    ComputeAnisotropicConstant(grid);

    //ComputePrandtlNumber(grid);

    double vistmax  = 0.0;
    double vistmin  = 1.0e30;
    double logsmall = log(SMALL);

    int nNegativeCr = 0;
    int nLowerBoundCr = 0;
    int nUpperBoundCr = 0;
    int nNegativeCi = 0;
    int nNegativePrt = 0;
    int nNegativeMu = 0;

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    int imin = 1;
    int jmin = 1;
    int kmin = 1;

    int nDim = GetDim();

    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int ist, ied, jst, jed, kst, ked;
        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble dudx = gradPrimtiveVarFaceXofLES(il, jl, kl, IU, iSurface);
                    RDouble dudy = gradPrimtiveVarFaceYofLES(il, jl, kl, IU, iSurface);
                    RDouble dudz = gradPrimtiveVarFaceZofLES(il, jl, kl, IU, iSurface);

                    RDouble dvdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IV, iSurface);

                    RDouble dwdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IW, iSurface);

                    RDouble dtdx = gradtx(il, jl, kl, ITT, iSurface);
                    RDouble dtdy = gradty(il, jl, kl, ITT, iSurface);
                    RDouble dtdz = gradtz(il, jl, kl, ITT, iSurface);

                    RDouble s11 = dudx;
                    RDouble s22 = dvdy;
                    RDouble s33 = dwdz;
                    RDouble s12 = half * (dudy + dvdx);
                    RDouble s13 = half * (dudz + dwdx);
                    RDouble s23 = half * (dvdz + dwdy);

                    RDouble sijsij = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                    RDouble ss = sqrt(two * sijsij);
                    RDouble skk = s11 + s22 + s33;

                    RDouble rho = 0.5 * (q(i, j, k, IR) + q(il, jl, kl, IR));
                    RDouble mul = 0.5 * (visl(i, j, k) + visl(il, jl, kl));
                    RDouble dampingFunction = wallFunction(il, jl, kl, iSurface);

                    RDouble cr = 0.5 * (anisotropicConstant(il, jl, kl) + anisotropicConstant(i, j, k));

                    if (cr < -0.01)
                    {
                        //cr = -0.01;
                        nLowerBoundCr = nLowerBoundCr + 1;
                    }

                    if (cr > 0.5)
                    {
                        //cr = 0.5;
                        nUpperBoundCr = nUpperBoundCr + 1;
                    }                   

                    RDouble mut = cr * dampingFunction * dampingFunction * rho * ss * refReNumber;

                    RDouble ci = 0.5 * (isotropicConstant(i, j, k) + isotropicConstant(il, jl, kl));
                    //ci = MIN(MAX(ci, 0.0), 0.01);

                    RDouble kSGS = 2.0 * ci * dampingFunction * dampingFunction * rho * ss * ss * refReNumber;

                    //oprt = 2.0 / (turbulentPrandtlNumber(il, jl, kl) + turbulentPrandtlNumber(i, j, k));
                    oprt = 1.111111;

                    RDouble hfxSGS = -mut * specificHeatAtConstantPressure * oprt * dtdx;
                    RDouble hfySGS = -mut * specificHeatAtConstantPressure * oprt * dtdy;
                    RDouble hfzSGS = -mut * specificHeatAtConstantPressure * oprt * dtdz;

                    subgridHeatFlux(il, jl, kl, 1, iSurface) = hfxSGS;
                    subgridHeatFlux(il, jl, kl, 2, iSurface) = hfySGS;
                    subgridHeatFlux(il, jl, kl, 3, iSurface) = hfzSGS;

                    if (mut < -mul)
                    {
                        nNegativeMu = nNegativeMu + 1;
                    }

                    if (turbViscousCutType == 1)
                    {
                        mut = MAX(mut - mul, 0.0);
                    }
                    else if (turbViscousCutType == 2)
                    {
                        mut = MAX(mut, 0.0);
                    }
                    else if (turbViscousCutType == 3)
                    {
                        mut = MAX(mut, -mul);
                    }

                    if (prtCutType == 1)
                    {
                        oprt = MAX(oprt, 0.0);
                    }

                    mut = MIN(eddyViscosityLimit, mut);

                    RDouble txxSGS = -2.0 * mut * s11 + 1.0/3.0 * kSGS;
                    RDouble tyySGS = -2.0 * mut * s22 + 1.0/3.0 * kSGS;
                    RDouble tzzSGS = -2.0 * mut * s33 + 1.0/3.0 * kSGS;
                    RDouble txySGS = -2.0 * mut * s12;
                    RDouble txzSGS = -2.0 * mut * s13;
                    RDouble tyzSGS = -2.0 * mut * s23;

                    subgridStress(il, jl, kl, 1, iSurface) = txxSGS;
                    subgridStress(il, jl, kl, 2, iSurface) = tyySGS;
                    subgridStress(il, jl, kl, 3, iSurface) = tzzSGS;
                    subgridStress(il, jl, kl, 4, iSurface) = txySGS;
                    subgridStress(il, jl, kl, 5, iSurface) = txzSGS;
                    subgridStress(il, jl, kl, 6, iSurface) = tyzSGS;     

                    if (vistmax < mut)
                    {
                        vistmax = mut;
                        imax = i;
                        jmax = j;
                        kmax = k;
                    }

                    if (vistmin > mut)
                    {
                        vistmin = mut;
                        imin = i;
                        jmin = j;
                        kmin = k;
                    }

                    if (cr < 0.0)
                    {
                        nNegativeCr = nNegativeCr + 1;
                    }

                    if (ci < 0.0)
                    {
                        nNegativeCi = nNegativeCi + 1;
                    }

                    if (oprt < 0.0)
                    {
                        nNegativePrt = nNegativePrt + 1;
                    }
                }
            }
        }

    }

    double vist_max = vistmax;
    double vist_min = vistmin;

    grid->UpdateData("vist_max", &vist_max, PHDOUBLE, 1);
    grid->UpdateData("vist_min", &vist_min, PHDOUBLE, 1);
    grid->UpdateData("nNegativeCr", &nNegativeCr, PHINT, 1);
    grid->UpdateData("nNegativeCi", &nNegativeCi, PHINT, 1);
    grid->UpdateData("nNegativePrt", &nNegativePrt, PHINT, 1);
    grid->UpdateData("nNegativeMu", &nNegativeMu, PHINT, 1);
    grid->UpdateData("nLowerBoundCr", &nLowerBoundCr, PHINT, 1);
    grid->UpdateData("nUpperBoundCr", &nUpperBoundCr, PHINT, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorVistmax)
        {
            cout << "vistmax = " << vistmax << " " << imax << " " << jmax << " " << kmax << "\n";
            cout << "vistmin = " << vistmin << " " << imin << " " << jmin << " " << kmin << "\n";
        }

        if (monitorNegativeConstant)
        {
            RDouble percentNegativeCr = nNegativeCr / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentNegativeMu = nNegativeMu / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentLowerBoundCr = nLowerBoundCr / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentUpperBoundCr = nUpperBoundCr / static_cast<RDouble> (ni * nj * nk * nDim);

            cout << "percentNegativeCr = " << percentNegativeCr << "\n";
            //cout << "nNegativeCi = " << nNegativeCi << "\n";
            //cout << "nNegativePrt = " << nNegativePrt << "\n";
            cout << "percentNegativeMu = " << percentNegativeMu << "\n";
            cout << "percentLowerBoundCr = " << percentLowerBoundCr << "\n";
            cout << "percentUpperBoundCr = " << percentUpperBoundCr << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void LESSolverStruct::DynamicSmagViscosityIncompressible(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_LESSolverStruct *parameters = GetControlParameters();

    //RDouble refReNumber  = parameters->GetRefReNumber();

    int monitorNegativeConstant = GlobalDataBase::GetIntParaFromDB("monitorNegativeConstant");

    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));

    ComputeAnisotropicConstant(grid);

    //ComputePrandtlNumber(grid);

    double vistmax  = 0.0;
    double vistmin  = 1.0e30;
    double logsmall = log(SMALL);

    int nNegativeCr = 0;
    int nLowerBoundCr = 0;
    int nUpperBoundCr = 0;

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble dampingFunction = wallFunction(i, j, k);

                RDouble cr = anisotropicConstant(i, j, k);

                if (cr < -0.01)
                {
                    //cr = -0.01;
                    nLowerBoundCr = nLowerBoundCr + 1;
                }

                if (cr > 0.5)
                {
                    //cr = 0.5;
                    nUpperBoundCr = nUpperBoundCr + 1;
                }                    

                if (cr < 0.0)
                {
                    nNegativeCr = nNegativeCr + 1;
                }

                anisotropicConstant(i, j, k) = cr * dampingFunction * dampingFunction;
            }
        }
    }

    RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();

    turbulentPrandtlNumber = prandtlTurbulence;

    isotropicConstant = 0.0;

    grid->UpdateData("nNegativeCr", &nNegativeCr, PHINT, 1);
    grid->UpdateData("nLowerBoundCr", &nLowerBoundCr, PHINT, 1);
    grid->UpdateData("nUpperBoundCr", &nUpperBoundCr, PHINT, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorNegativeConstant)
        {
            RDouble percentNegativeCr = nNegativeCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentLowerBoundCr = nLowerBoundCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentUpperBoundCr = nUpperBoundCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));

            cout << "percentNegativeCr = " << percentNegativeCr << "\n";
            cout << "percentLowerBoundCr = " << percentLowerBoundCr << "\n";
            cout << "percentUpperBoundCr = " << percentUpperBoundCr << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void LESSolverStruct::ComputeTestFilteredQ(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    testFilteredQ(I, J, K, Range(IR, IW)) = q(I, J, K, Range(IR, IW));

    //! filter Q in I direction
    if (filterDirection[0] == 1) 
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredQ(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredQ(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredQ(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter Q in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredQ(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredQ(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredQ(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter Q in K direction
    if (filterDirection[2] == 1) 
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredQ(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredQ(i, j, K, n);
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredQ(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredQfd(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    testFilteredQ(I, J, K, Range(IR, IW)) = q(I, J, K, Range(IR, IW));

    testFilteredQ(I, J, K, Range(IU, IW)) = q(I, J, K, IR) * q(I, J, K, Range(IU, IW));

    //! filter Q in I direction
    if (filterDirection[0] == 1) 
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredQ(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredQ(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredQ(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter Q in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredQ(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredQ(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredQ(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter Q in K direction
    if (filterDirection[2] == 1) 
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredQ(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredQ(i, j, K, n);
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredQ(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    testFilteredQ(I, J, K, Range(IU, IW)) = testFilteredQ(I, J, K, Range(IU, IW)) / testFilteredQ(I, J, K, IR);

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

/*
void LESSolverStruct::ComputeTestFilteredDensity(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D &testFilteredDensity = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredDensity"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    testFilteredDensity(I, J, K) = q(I, J, K, IR);

    //! filter density in I direction
    if (filterDirection[0] == 1) 
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                (*testFilteredVariableX)(I) = testFilteredDensity(I, j, k);
                (*filteredVariableX)(I) = testFilteredDensity(I, j, k);
                TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                testFilteredDensity(I, j, k) = (*testFilteredVariableX)(I);
            }
        }
    }

    //! filter density in J direction
    if (filterDirection[1] == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*testFilteredVariableY)(J) = testFilteredDensity(i, J, k);
                (*filteredVariableY)(J) = testFilteredDensity(i, J, k);
                TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                testFilteredDensity(i, J, k) = (*testFilteredVariableY)(J);
            }
        }
    }

    //! filter density in K direction
    if (filterDirection[2] == 1)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*testFilteredVariableZ)(K) = testFilteredDensity(i, j, K);
                (*filteredVariableZ)(K) = testFilteredDensity(i, j, K);
                TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                testFilteredDensity(i, j, K) = (*testFilteredVariableZ)(K);
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}
*/

void LESSolverStruct::ComputeTestFilteredTemperature(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &temperatures = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("t"));

    RDouble3D &testFilteredTemperature = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredTemperature"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    testFilteredTemperature(I, J, K) = temperatures(I, J, K, ITT);

    //! filter density in I direction
    if (filterDirection[0] == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                (*testFilteredVariableX)(I) = testFilteredTemperature(I, j, k);
                (*filteredVariableX)(I) = testFilteredTemperature(I, j, k);
                TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                testFilteredTemperature(I, j, k) = (*testFilteredVariableX)(I);
            }
        }
    }

    //! filter density in J direction
    if (filterDirection[1] == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*testFilteredVariableY)(J) = testFilteredTemperature(i, J, k);
                (*filteredVariableY)(J) = testFilteredTemperature(i, J, k);
                TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                testFilteredTemperature(i, J, k) = (*testFilteredVariableY)(J);
            }
        }
    }

    //! filter density in K direction
    if (filterDirection[2] == 1)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*testFilteredVariableZ)(K) = testFilteredTemperature(i, j, K);
                (*filteredVariableZ)(K) = testFilteredTemperature(i, j, K);
                TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                testFilteredTemperature(i, j, K) = (*testFilteredVariableZ)(K);
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeStrainRateTensor(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));
    //RDouble3D &strainRateMagnitude = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("strainRateMagnitude"));

    /*
    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    strainRateTensor(I, J, K, 1) = gradUVWTCellCenterX(I, J, K, 0);
    strainRateTensor(I, J, K, 2) = gradUVWTCellCenterY(I, J, K, 1);
    strainRateTensor(I, J, K, 3) = gradUVWTCellCenterZ(I, J, K, 2);
    strainRateTensor(I, J, K, 4) = 0.5 * (gradUVWTCellCenterY(I, J, K, 0) + gradUVWTCellCenterX(I, J, K, 1));
    strainRateTensor(I, J, K, 5) = 0.5 * (gradUVWTCellCenterZ(I, J, K, 0) + gradUVWTCellCenterX(I, J, K, 2));
    strainRateTensor(I, J, K, 6) = 0.5 * (gradUVWTCellCenterY(I, J, K, 2) + gradUVWTCellCenterZ(I, J, K, 1));

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble s11 = strainRateTensor(i, j, k, 1);
                RDouble s22 = strainRateTensor(i, j, k, 2);
                RDouble s33 = strainRateTensor(i, j, k, 3);
                RDouble s12 = strainRateTensor(i, j, k, 4);
                RDouble s13 = strainRateTensor(i, j, k, 5);
                RDouble s23 = strainRateTensor(i, j, k, 6);

                RDouble s = 2.0 * (s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23));
                strainRateMagnitude(i, j, k) = sqrt(s);
            }
        }
    }
    */

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D gradux(I, J, K, fortranArray);
    RDouble3D graduy(I, J, K, fortranArray);
    RDouble3D graduz(I, J, K, fortranArray);

    RDouble3D gradvx(I, J, K, fortranArray);
    RDouble3D gradvy(I, J, K, fortranArray);
    RDouble3D gradvz(I, J, K, fortranArray);

    RDouble3D gradwx(I, J, K, fortranArray);
    RDouble3D gradwy(I, J, K, fortranArray);
    RDouble3D gradwz(I, J, K, fortranArray);

    GradCenter(grid, q, gradux, graduy, graduz, IU);
    GradCenter(grid, q, gradvx, gradvy, gradvz, IV);
    GradCenter(grid, q, gradwx, gradwy, gradwz, IW);

    strainRateTensor(I, J, K, 1) = gradux(I, J, K);
    strainRateTensor(I, J, K, 2) = gradvy(I, J, K);
    strainRateTensor(I, J, K, 3) = gradwz(I, J, K);
    strainRateTensor(I, J, K, 4) = 0.5 * (graduy(I, J, K) + gradvx(I, J, K));
    strainRateTensor(I, J, K, 5) = 0.5 * (graduz(I, J, K) + gradwx(I, J, K));
    strainRateTensor(I, J, K, 6) = 0.5 * (gradwy(I, J, K) + gradvz(I, J, K));

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble s11 = strainRateTensor(i, j, k, 1);
                RDouble s22 = strainRateTensor(i, j, k, 2);
                RDouble s33 = strainRateTensor(i, j, k, 3);
                RDouble s12 = strainRateTensor(i, j, k, 4);
                RDouble s13 = strainRateTensor(i, j, k, 5);
                RDouble s23 = strainRateTensor(i, j, k, 6);

                RDouble s = 2.0 * (s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23));
                strainRateTensor(i, j, k, 0) = sqrt(s);
            }
        }
    }

    GhostCell3D(strainRateTensor, ni, nj, nk, 7);
}

void LESSolverStruct::ComputeStrainRateTensorforStructHighOrder(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));
    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D gradux(I, J, K, fortranArray);
    RDouble3D graduy(I, J, K, fortranArray);
    RDouble3D graduz(I, J, K, fortranArray);

    RDouble3D gradvx(I, J, K, fortranArray);
    RDouble3D gradvy(I, J, K, fortranArray);
    RDouble3D gradvz(I, J, K, fortranArray);

    RDouble3D gradwx(I, J, K, fortranArray);
    RDouble3D gradwy(I, J, K, fortranArray);
    RDouble3D gradwz(I, J, K, fortranArray);

    GradCenter(grid, q, gradux, graduy, graduz, IU);
    GradCenter(grid, q, gradvx, gradvy, gradvz, IV);
    GradCenter(grid, q, gradwx, gradwy, gradwz, IW);

    strainRateTensor(I, J, K, 1) = gradux(I, J, K);
    strainRateTensor(I, J, K, 2) = gradvy(I, J, K);
    strainRateTensor(I, J, K, 3) = gradwz(I, J, K);
    strainRateTensor(I, J, K, 4) = 0.5 * (graduy(I, J, K) + gradvx(I, J, K));
    strainRateTensor(I, J, K, 5) = 0.5 * (graduz(I, J, K) + gradwx(I, J, K));
    strainRateTensor(I, J, K, 6) = 0.5 * (gradwy(I, J, K) + gradvz(I, J, K));

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble s11 = strainRateTensor(i, j, k, 1);
                RDouble s22 = strainRateTensor(i, j, k, 2);
                RDouble s33 = strainRateTensor(i, j, k, 3);
                RDouble s12 = strainRateTensor(i, j, k, 4);
                RDouble s13 = strainRateTensor(i, j, k, 5);
                RDouble s23 = strainRateTensor(i, j, k, 6);

                RDouble s = 2.0 * (s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23));
                strainRateTensor(i, j, k, 0) = sqrt(s);
            }
        }
    }

    GhostCell3D(strainRateTensor, ni, nj, nk, 7);
}

void LESSolverStruct::ComputeGradT(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    RDouble4D &temperatures = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("t"));
    RDouble4D &gradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradT"));

    RDouble3D gradTx(I, J, K, fortranArray);
    RDouble3D gradTy(I, J, K, fortranArray);
    RDouble3D gradTz(I, J, K, fortranArray);

    GradCenter(grid, temperatures, gradTx, gradTy, gradTz, ITT);

    GhostCell3D(gradTx, ni, nj, nk);
    GhostCell3D(gradTy, ni, nj, nk);
    GhostCell3D(gradTz, ni, nj, nk);

    gradT(I, J, K, 0) = gradTx(I, J, K);
    gradT(I, J, K, 1) = gradTy(I, J, K);
    gradT(I, J, K, 2) = gradTz(I, J, K);
}

void LESSolverStruct::ComputeTestFilteredStrainRateTensor(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range N(0, 6);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));
    RDouble4D &testFilteredStrainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredStrainRateTensor"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    testFilteredStrainRateTensor(I, J, K, N) = strainRateTensor(I, J, K, N);

    //! filter Sij in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredStrainRateTensor(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredStrainRateTensor(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredStrainRateTensor(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    GhostCell3D(testFilteredStrainRateTensor, ni, nj, nk, 7);

    //! filter Sij in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredStrainRateTensor(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredStrainRateTensor(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredStrainRateTensor(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    GhostCell3D(testFilteredStrainRateTensor, ni, nj, nk, 7);

    //! filter Sij in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredStrainRateTensor(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredStrainRateTensor(i, j, K, n);
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredStrainRateTensor(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    GhostCell3D(testFilteredStrainRateTensor, ni, nj, nk, 7);

    //! Compute |S|
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble tFS11 = testFilteredStrainRateTensor(i, j, k, 1);
                RDouble tFS22 = testFilteredStrainRateTensor(i, j, k, 2);
                RDouble tFS33 = testFilteredStrainRateTensor(i, j, k, 3);
                RDouble tFS12 = testFilteredStrainRateTensor(i, j, k, 4);
                RDouble tFS13 = testFilteredStrainRateTensor(i, j, k, 5);
                RDouble tFS23 = testFilteredStrainRateTensor(i, j, k, 6);

                RDouble tFS = 2.0 * (tFS11 * tFS11 + tFS22 * tFS22 + tFS33 * tFS33 + 2.0 * (tFS12 * tFS12 + tFS13 * tFS13 + tFS23 * tFS23));
                testFilteredStrainRateTensor(i, j, k, 0) = sqrt(tFS);
            }
        }
    }

    GhostCell3D(testFilteredStrainRateTensor, ni, nj, nk, 1);

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredRhoUiUj(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();    
    int *filterDirection = parameters->GetFilterDirection();
    int subgridScaleModel  = parameters->GetSubgridScaleModel();

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &testFilteredRhoUiUj = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUiUj"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    if (subgridScaleModel == 2)
    {
        testFilteredRhoUiUj(I, J, K, 1) = q(I, J, K, IR) * q(I, J, K, IU) * q(I, J, K, IU);
        testFilteredRhoUiUj(I, J, K, 2) = q(I, J, K, IR) * q(I, J, K, IV) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 3) = q(I, J, K, IR) * q(I, J, K, IW) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 4) = q(I, J, K, IR) * q(I, J, K, IU) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 5) = q(I, J, K, IR) * q(I, J, K, IU) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 6) = q(I, J, K, IR) * q(I, J, K, IV) * q(I, J, K, IW);
    }
    else if (subgridScaleModel == 4)
    {
        testFilteredRhoUiUj(I, J, K, 1) = q(I, J, K, IU) * q(I, J, K, IU);
        testFilteredRhoUiUj(I, J, K, 2) = q(I, J, K, IV) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 3) = q(I, J, K, IW) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 4) = q(I, J, K, IU) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 5) = q(I, J, K, IU) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 6) = q(I, J, K, IV) * q(I, J, K, IW);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }

    //! filter RhoUiUj in I direction
    //testFilteredRhoUiUj(ist-1, Range(jst, jed), Range(kst, ked), Range(1, 6)) = testFilteredRhoUiUj(ist, Range(jst, jed), Range(kst, ked), Range(1, 6));
    //testFilteredRhoUiUj(ied+1, Range(jst, jed), Range(kst, ked), Range(1, 6)) = testFilteredRhoUiUj(ied, Range(jst, jed), Range(kst, ked), Range(1, 6));

    RDouble tFRhoUU = testFilteredRhoUiUj(1, 1, 1, 1);
    RDouble tFRhoVV = testFilteredRhoUiUj(1, 1, 1, 2);
    RDouble tFRhoWW = testFilteredRhoUiUj(1, 1, 1, 3);
    RDouble tFRhoUV = testFilteredRhoUiUj(1, 1, 1, 4);
    RDouble tFRhoUW = testFilteredRhoUiUj(1, 1, 1, 5);
    RDouble tFRhoVW = testFilteredRhoUiUj(1, 1, 1, 6);

    if (filterDirection[0] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredRhoUiUj(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredRhoUiUj(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredRhoUiUj(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    tFRhoUU = testFilteredRhoUiUj(1, 1, 1, 1);
    tFRhoVV = testFilteredRhoUiUj(1, 1, 1, 2);
    tFRhoWW = testFilteredRhoUiUj(1, 1, 1, 3);
    tFRhoUV = testFilteredRhoUiUj(1, 1, 1, 4);
    tFRhoUW = testFilteredRhoUiUj(1, 1, 1, 5);
    tFRhoVW = testFilteredRhoUiUj(1, 1, 1, 6);

    //! filter RhoUiUj in J direction
    //testFilteredRhoUiUj(Range(ist, ied), jst-1, Range(kst, ked), Range(1, 6)) = testFilteredRhoUiUj(Range(ist, ied), jst, Range(kst, ked), Range(1, 6));
    //testFilteredRhoUiUj(Range(ist, ied), jed+1, Range(kst, ked), Range(1, 6)) = testFilteredRhoUiUj(Range(ist, ied), jed, Range(kst, ked), Range(1, 6));

    if (filterDirection[1] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredRhoUiUj(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredRhoUiUj(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredRhoUiUj(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    tFRhoUU = testFilteredRhoUiUj(1, 1, 1, 1);
    tFRhoVV = testFilteredRhoUiUj(1, 1, 1, 2);
    tFRhoWW = testFilteredRhoUiUj(1, 1, 1, 3);
    tFRhoUV = testFilteredRhoUiUj(1, 1, 1, 4);
    tFRhoUW = testFilteredRhoUiUj(1, 1, 1, 5);
    tFRhoVW = testFilteredRhoUiUj(1, 1, 1, 6);

    //! filter RhoUiUj in K direction
    //testFilteredRhoUiUj(Range(ist, ied), Range(jst, jed), kst-1, Range(1, 6)) = testFilteredRhoUiUj(Range(ist, ied), Range(jst, jed), kst, Range(1, 6));
    //testFilteredRhoUiUj(Range(ist, ied), Range(jst, jed), ked+1, Range(1, 6)) = testFilteredRhoUiUj(Range(ist, ied), Range(jst, jed), ked, Range(1, 6));

    if (filterDirection[2] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredRhoUiUj(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredRhoUiUj(i, j, K, n);                
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredRhoUiUj(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeLeonardStress(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    RDouble4D &leonardStress = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardStress"));
    RDouble4D &testFilteredRhoUiUj = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUiUj"));
    //RDouble3D &testFilteredDensity = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredDensity"));
    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));

    ComputeTestFilteredQ(grid);

    Param_LESSolverStruct * parameters = GetControlParameters();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    if (subgridScaleModel == 4)
    {
        testFilteredQ(I, J, K, 0) = 1.0;
    }

    ComputeTestFilteredRhoUiUj(grid);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble tFRho = testFilteredQ(i, j, k, 0);

                RDouble tFU = testFilteredQ(i, j, k, 1);
                RDouble tFV = testFilteredQ(i, j, k, 2);
                RDouble tFW = testFilteredQ(i, j, k, 3);

                RDouble tFRhoUU = testFilteredRhoUiUj(i, j, k, 1);
                RDouble tFRhoVV = testFilteredRhoUiUj(i, j, k, 2);
                RDouble tFRhoWW = testFilteredRhoUiUj(i, j, k, 3);
                RDouble tFRhoUV = testFilteredRhoUiUj(i, j, k, 4);
                RDouble tFRhoUW = testFilteredRhoUiUj(i, j, k, 5);
                RDouble tFRhoVW = testFilteredRhoUiUj(i, j, k, 6);

                leonardStress(i, j, k, 1) = tFRhoUU - tFRho * tFU * tFU;
                leonardStress(i, j, k, 2) = tFRhoVV - tFRho * tFV * tFV;
                leonardStress(i, j, k, 3) = tFRhoWW - tFRho * tFW * tFW;
                leonardStress(i, j, k, 4) = tFRhoUV - tFRho * tFU * tFV;
                leonardStress(i, j, k, 5) = tFRhoUW - tFRho * tFU * tFW;
                leonardStress(i, j, k, 6) = tFRhoVW - tFRho * tFV * tFW;
            }
        }
    }
}

void LESSolverStruct::ComputeModeledStress(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range N(1, 6);

    RDouble4D &modeledStress = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledStress"));
    RDouble4D &testFilteredAlphaIJ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredAlphaIJ"));
    RDouble4D &betaIJ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("betaIJ"));

    ComputeTestFilteredAlphaIJ(grid);

    ComputeBetaIJ(grid);

    modeledStress(I, J, K, N) = betaIJ(I, J, K, N) - testFilteredAlphaIJ(I, J, K, N); 
}

void LESSolverStruct::ComputeTestFilteredAlphaIJ(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));

    RDouble4D &testFilteredAlphaIJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("testFilteredAlphaIJ"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    if (subgridScaleModel == 2)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble rho = q(i, j, k, IR);
                    RDouble s = strainRateTensor(i, j, k, 0);
                    RDouble s11 = strainRateTensor(i, j, k, 1);
                    RDouble s22 = strainRateTensor(i, j, k, 2);
                    RDouble s33 = strainRateTensor(i, j, k, 3);
                    RDouble s12 = strainRateTensor(i, j, k, 4);
                    RDouble s13 = strainRateTensor(i, j, k, 5);
                    RDouble s23 = strainRateTensor(i, j, k, 6);

                    testFilteredAlphaIJ(i, j, k, 1) = -2.0 * rho * s * (s11 - 1.0/3.0 * (s11 + s22 + s33));
                    testFilteredAlphaIJ(i, j, k, 2) = -2.0 * rho * s * (s22 - 1.0/3.0 * (s11 + s22 + s33));
                    testFilteredAlphaIJ(i, j, k, 3) = -2.0 * rho * s * (s33 - 1.0/3.0 * (s11 + s22 + s33));
                    testFilteredAlphaIJ(i, j, k, 4) = -2.0 * rho * s * s12;
                    testFilteredAlphaIJ(i, j, k, 5) = -2.0 * rho * s * s13;
                    testFilteredAlphaIJ(i, j, k, 6) = -2.0 * rho * s * s23;

                    testFilteredAlphaIJ(i, j, k, 0) = 2.0 * rho * s * s;
                }
            }
        }
    }
    else if (subgridScaleModel == 4)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble s = strainRateTensor(i, j, k, 0);
                    RDouble s11 = strainRateTensor(i, j, k, 1);
                    RDouble s22 = strainRateTensor(i, j, k, 2);
                    RDouble s33 = strainRateTensor(i, j, k, 3);
                    RDouble s12 = strainRateTensor(i, j, k, 4);
                    RDouble s13 = strainRateTensor(i, j, k, 5);
                    RDouble s23 = strainRateTensor(i, j, k, 6);

                    testFilteredAlphaIJ(i, j, k, 1) = -2.0 * s * s11;
                    testFilteredAlphaIJ(i, j, k, 2) = -2.0 * s * s22;
                    testFilteredAlphaIJ(i, j, k, 3) = -2.0 * s * s33;
                    testFilteredAlphaIJ(i, j, k, 4) = -2.0 * s * s12;
                    testFilteredAlphaIJ(i, j, k, 5) = -2.0 * s * s13;
                    testFilteredAlphaIJ(i, j, k, 6) = -2.0 * s * s23;

                    testFilteredAlphaIJ(i, j, k, 0) = 2.0 * s * s;
                }
            }
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }

    GhostCell3D(testFilteredAlphaIJ, ni, nj, nk, 7);

    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    //! filter AlphaIJ in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 0; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredAlphaIJ(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredAlphaIJ(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredAlphaIJ(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    GhostCell3D(testFilteredAlphaIJ, ni, nj, nk, 7);

    //! filter AlphaIJ in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 0; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredAlphaIJ(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredAlphaIJ(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredAlphaIJ(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    GhostCell3D(testFilteredAlphaIJ, ni, nj, nk, 7);

    //! filter AlphaIJ in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 0; n <= 6; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredAlphaIJ(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredAlphaIJ(i, j, K, n);                
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredAlphaIJ(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    GhostCell3D(testFilteredAlphaIJ, ni, nj, nk, 7);

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeBetaIJ(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();
    RDouble testFilterScale = parameters->GetTestFilterScale();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    //RDouble3D &testFilteredDensity = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredDensity"));
    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));
    RDouble4D &testFilteredStrainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredStrainRateTensor"));

    RDouble4D &betaIJ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("betaIJ"));

    ComputeTestFilteredStrainRateTensor(grid);

    if (subgridScaleModel == 2)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble rho = testFilteredQ(i, j, k, 0);
                    RDouble s = testFilteredStrainRateTensor(i, j, k, 0);
                    RDouble s11 = testFilteredStrainRateTensor(i, j, k, 1);
                    RDouble s22 = testFilteredStrainRateTensor(i, j, k, 2);
                    RDouble s33 = testFilteredStrainRateTensor(i, j, k, 3);
                    RDouble s12 = testFilteredStrainRateTensor(i, j, k, 4);
                    RDouble s13 = testFilteredStrainRateTensor(i, j, k, 5);
                    RDouble s23 = testFilteredStrainRateTensor(i, j, k, 6);

                    betaIJ(i, j, k, 1) = -2.0 * testFilterScale * testFilterScale * rho * s * (s11 - 1.0/3.0 * (s11 + s22 + s33));
                    betaIJ(i, j, k, 2) = -2.0 * testFilterScale * testFilterScale * rho * s * (s22 - 1.0/3.0 * (s11 + s22 + s33));
                    betaIJ(i, j, k, 3) = -2.0 * testFilterScale * testFilterScale * rho * s * (s33 - 1.0/3.0 * (s11 + s22 + s33));
                    betaIJ(i, j, k, 4) = -2.0 * testFilterScale * testFilterScale * rho * s * s12;
                    betaIJ(i, j, k, 5) = -2.0 * testFilterScale * testFilterScale * rho * s * s13;
                    betaIJ(i, j, k, 6) = -2.0 * testFilterScale * testFilterScale * rho * s * s23;
                }
            }
        }
    }
    else if (subgridScaleModel == 4)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble s = testFilteredStrainRateTensor(i, j, k, 0);
                    RDouble s11 = testFilteredStrainRateTensor(i, j, k, 1);
                    RDouble s22 = testFilteredStrainRateTensor(i, j, k, 2);
                    RDouble s33 = testFilteredStrainRateTensor(i, j, k, 3);
                    RDouble s12 = testFilteredStrainRateTensor(i, j, k, 4);
                    RDouble s13 = testFilteredStrainRateTensor(i, j, k, 5);
                    RDouble s23 = testFilteredStrainRateTensor(i, j, k, 6);

                    betaIJ(i, j, k, 1) = -2.0 * testFilterScale * testFilterScale * s * s11;
                    betaIJ(i, j, k, 2) = -2.0 * testFilterScale * testFilterScale * s * s22;
                    betaIJ(i, j, k, 3) = -2.0 * testFilterScale * testFilterScale * s * s33;
                    betaIJ(i, j, k, 4) = -2.0 * testFilterScale * testFilterScale * s * s12;
                    betaIJ(i, j, k, 5) = -2.0 * testFilterScale * testFilterScale * s * s13;
                    betaIJ(i, j, k, 6) = -2.0 * testFilterScale * testFilterScale * s * s23;
                }
            }
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }
}

void LESSolverStruct::ComputeAnisotropicConstant(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble4D &leonardStress = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardStress"));
    RDouble4D &modeledStress = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledStress"));

    ComputeLeonardStress(grid);

    ComputeModeledStress(grid);

    RDouble3D *LijMij = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *MijMij = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *LmmMnn = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *LijMijAvg = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *MijMijAvg = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *LmmMnnAvg = new RDouble3D(I, J, K, fortranArray);

    if (subgridScaleModel == 2)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble l11 = leonardStress(i, j, k, 1);
                    RDouble l22 = leonardStress(i, j, k, 2);
                    RDouble l33 = leonardStress(i, j, k, 3);
                    RDouble l12 = leonardStress(i, j, k, 4);
                    RDouble l13 = leonardStress(i, j, k, 5);
                    RDouble l23 = leonardStress(i, j, k, 6);

                    RDouble m11 = modeledStress(i, j, k, 1);
                    RDouble m22 = modeledStress(i, j, k, 2);
                    RDouble m33 = modeledStress(i, j, k, 3);
                    RDouble m12 = modeledStress(i, j, k, 4);
                    RDouble m13 = modeledStress(i, j, k, 5);
                    RDouble m23 = modeledStress(i, j, k, 6);

                    RDouble lm = l11 * m11 + l22 * m22 + l33 * m33 + 2.0 * (l12 * m12 + l13 * m13 + l23 * m23);
                    RDouble mm = m11 * m11 + m22 * m22 + m33 * m33 + 2.0 * (m12 * m12 + m13 * m13 + m23 * m23);
                    RDouble lmmn = (l11 + l22 + l33) * (m11 + m22 + m33);

                    (*LijMij)(i, j, k) = lm;
                    (*MijMij)(i, j, k) = mm;
                    (*LmmMnn)(i, j, k) = lmmn;
                }
            }
        }
    }
    else if (subgridScaleModel == 4)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble l11 = leonardStress(i, j, k, 1);
                    RDouble l22 = leonardStress(i, j, k, 2);
                    RDouble l33 = leonardStress(i, j, k, 3);
                    RDouble l12 = leonardStress(i, j, k, 4);
                    RDouble l13 = leonardStress(i, j, k, 5);
                    RDouble l23 = leonardStress(i, j, k, 6);

                    RDouble m11 = modeledStress(i, j, k, 1);
                    RDouble m22 = modeledStress(i, j, k, 2);
                    RDouble m33 = modeledStress(i, j, k, 3);
                    RDouble m12 = modeledStress(i, j, k, 4);
                    RDouble m13 = modeledStress(i, j, k, 5);
                    RDouble m23 = modeledStress(i, j, k, 6);

                    RDouble lm = l11 * m11 + l22 * m22 + l33 * m33 + 2.0 * (l12 * m12 + l13 * m13 + l23 * m23);
                    RDouble mm = m11 * m11 + m22 * m22 + m33 * m33 + 2.0 * (m12 * m12 + m13 * m13 + m23 * m23);
                    RDouble lmmn = (l11 + l22 + l33) * (m11 + m22 + m33);

                    (*LijMij)(i, j, k) = lm;
                    (*MijMij)(i, j, k) = mm;
                    (*LmmMnn)(i, j, k) = 0.0;
                }
            }
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }

    (*LijMijAvg)(I, J, K) = (*LijMij)(I, J, K);
    (*MijMijAvg)(I, J, K) = (*MijMij)(I, J, K);
    (*LmmMnnAvg)(I, J, K) = (*LmmMnn)(I, J, K);

    int *averageDirection = parameters->GetAverageDirection();
    
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    if (averageDirection[0] == 0 && averageDirection[1] == 0 && averageDirection[2] == 0)
    {
        int averageWidth = parameters->GetAverageWidth();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble LMAvg = 0.0;
                    RDouble MMAvg = 0.0;
                    RDouble LmMnAvg = 0.0;

                    int nAvg = 0;

                    for (int ka = -averageWidth; ka <= averageWidth; ++ ka)
                    {
                        for (int ja = -averageWidth; ja <= averageWidth; ++ ja)
                        {
                            for (int ia = -averageWidth; ia <= averageWidth; ++ ia)
                            {
                                int kk = k + ka;
                                int jj = j + ja;
                                int ii = i + ia;

                                bool ifAverage = (ii >= ist && ii <= ied) && (jj >= jst && jj <= jed) && (kk >= kst && kk <= ked);

                                if (ifAverage)
                                {
                                    LMAvg = LMAvg + (*LijMij)(ii, jj, kk);
                                    MMAvg = MMAvg + (*MijMij)(ii, jj, kk);
                                    LmMnAvg = LmMnAvg + (*LmmMnn)(ii, jj, kk);
                                    nAvg = nAvg + 1;
                                }
                            }
                        }
                    }

                    (*LijMijAvg)(i, j, k) = LMAvg / static_cast<RDouble>(nAvg);
                    (*MijMijAvg)(i, j, k) = MMAvg / static_cast<RDouble>(nAvg);
                    (*LmmMnnAvg)(i, j, k) = LmMnAvg / static_cast<RDouble>(nAvg);
                }
            }
        }
    }
    else
    {
        //! average in I direction
        if (averageDirection[0] != 0)
        {
            int nAvg = ied - ist + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    RDouble LMAvg = 0.0;
                    RDouble MMAvg = 0.0;
                    RDouble LmMnAvg = 0.0;

                    for (int i = ist; i <= ied; ++ i)
                    {
                        LMAvg = LMAvg + (*LijMij)(i, j, k);
                        MMAvg = MMAvg + (*MijMij)(i, j, k);
                        LmMnAvg = LmMnAvg + (*LmmMnn)(i, j, k);
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*LijMijAvg)(i, j, k) = LMAvg / static_cast <RDouble>(nAvg);
                        (*MijMijAvg)(i, j, k) = MMAvg / static_cast <RDouble>(nAvg);
                        (*LmmMnnAvg)(i, j, k) = LmMnAvg / static_cast <RDouble>(nAvg);
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*LijMij)(i, j, k) = (*LijMijAvg)(i, j, k);
                        (*MijMij)(i, j, k) = (*MijMijAvg)(i, j, k);
                        (*LmmMnn)(i, j, k) = (*LmmMnnAvg)(i, j, k);
                    }
                }
            }
        }

        RDouble lm = (*LijMijAvg)(1, 1, 1);
        RDouble mm = (*MijMijAvg)(1, 1, 1);

        //! average in J direction
        if (averageDirection[1] != 0)
        {
            int nAvg = jed - jst + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble LMAvg = 0.0;
                    RDouble MMAvg = 0.0;
                    RDouble LmMnAvg = 0.0;

                    for (int j = jst; j <= jed; ++ j)
                    {
                        RDouble a = (*LijMij)(i, j, k);
                        RDouble b = (*MijMij)(i, j, k);
                        LMAvg = LMAvg + (*LijMij)(i, j, k);
                        MMAvg = MMAvg + (*MijMij)(i, j, k);
                        LmMnAvg = LmMnAvg + (*LmmMnn)(i, j, k);
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*LijMijAvg)(i, j, k) = LMAvg / static_cast <RDouble>(nAvg);
                        (*MijMijAvg)(i, j, k) = MMAvg / static_cast <RDouble>(nAvg);
                        (*LmmMnnAvg)(i, j, k) = LmMnAvg / static_cast <RDouble>(nAvg);
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*LijMij)(i, j, k) = (*LijMijAvg)(i, j, k);
                        (*MijMij)(i, j, k) = (*MijMijAvg)(i, j, k);
                        (*LmmMnn)(i, j, k) = (*LmmMnnAvg)(i, j, k);
                    }
                }
            }
        }

        lm = (*LijMijAvg)(1, 1, 1);
        mm = (*MijMijAvg)(1, 1, 1);

        //! average in K direction
        if (averageDirection[2] != 0)
        {
            int nAvg = ked - kst + 1;

            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble LMAvg = 0.0;
                    RDouble MMAvg = 0.0;
                    RDouble LmMnAvg = 0.0;

                    for (int k = kst; k <= ked; ++ k)
                    {
                        LMAvg = LMAvg + (*LijMij)(i, j, k);
                        MMAvg = MMAvg + (*MijMij)(i, j, k);
                        LmMnAvg = LmMnAvg + (*LmmMnn)(i, j, k);
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*LijMijAvg)(i, j, k) = LMAvg / static_cast <RDouble>(nAvg);
                        (*MijMijAvg)(i, j, k) = MMAvg / static_cast <RDouble>(nAvg);
                        (*LmmMnnAvg)(i, j, k) = LmMnAvg / static_cast <RDouble>(nAvg);
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*LijMij)(i, j, k) = (*LijMijAvg)(i, j, k);
                        (*MijMij)(i, j, k) = (*MijMijAvg)(i, j, k);
                        (*LmmMnn)(i, j, k) = (*LmmMnnAvg)(i, j, k);
                    }
                }
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble lm = (*LijMijAvg)(i, j, k);
                RDouble mm = (*MijMijAvg)(i, j, k);
                RDouble lmMn = (*LmmMnnAvg)(i, j, k);

                anisotropicConstant(i, j, k) = (lm - 1.0/3.0 * lmMn) / mm; 
            }
        }
    }

    GhostCell3D(anisotropicConstant, ni, nj, nk);

    delete LijMij;
    delete MijMij;
    delete LmmMnn;
    delete LijMijAvg;
    delete MijMijAvg;
    delete LmmMnnAvg;
}

void LESSolverStruct::ComputeIsotropicConstant(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    RDouble testFilterScale = parameters->GetTestFilterScale();
    int *averageDirection = parameters->GetAverageDirection();

    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));
    RDouble4D &leonardStress = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardStress"));
    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));
    RDouble4D &testFilteredAlphaIJ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredAlphaIJ"));
    RDouble4D &testFilteredStrainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredStrainRateTensor"));

    RDouble3D *numerator = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *denominator = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *numeratorAvg = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *denominatorAvg = new RDouble3D(I, J, K, fortranArray);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble l11 = leonardStress(i, j, k, 1);
                RDouble l22 = leonardStress(i, j, k, 2);
                RDouble l33 = leonardStress(i, j, k, 3);
                (*numerator)(i, j, k) = l11 + l22 + l33;

                RDouble tFRho = testFilteredQ(i, j, k, 0);
                RDouble testFilteredStranRateMagnitude = testFilteredStrainRateTensor(i, j, k, 0);
                RDouble testFilteredAlpha = testFilteredAlphaIJ(i, j, k, 0);
                (*denominator)(i, j, k) = 2.0 * testFilterScale * testFilterScale * tFRho * testFilteredStranRateMagnitude * testFilteredStranRateMagnitude - testFilteredAlpha;
            }
        }
    }

    //(*numerator)(I, J, K) = leonardStress(I, J, K, 1) + leonardStress(I, J, K, 2) + leonardStress(I, J, K, 3);
    //(*denominator)(I, J, K) = 2.0 * testFilterScale * testFilterScale * testFilteredQ(I, J, K, 0) * testFilteredStrainRateTensor(I, J, K, 0) * testFilteredStrainRateTensor(I, J, K, 0) - testFilteredAlphaIJ(I, J, K, 0);

    (*numeratorAvg)(I, J, K) = (*numerator)(I, J, K);
    (*denominatorAvg)(I, J, K) = (*denominator)(I, J, K);

    if (averageDirection[0] == 0 && averageDirection[1] == 0 && averageDirection[2] == 0)
    {
        int averageWidth = parameters->GetAverageWidth();
        
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    int nAvg = 0;

                    for (int ka = -averageWidth; ka <= averageWidth; ++ ka)
                    {
                        for (int ja = -averageWidth; ja <= averageWidth; ++ ja)
                        {
                            for (int ia = -averageWidth; ia <= averageWidth; ++ ia)
                            {
                                int kk = k + ka;
                                int jj = j + ja;
                                int ii = i + ia;

                                bool ifAverage = (ii >= ist && ii <= ied) && (jj >= jst && jj <= jed) && (kk >= kst && kk <= ked);

                                if (ifAverage)
                                {
                                    num = num + (*numerator)(ii, jj, kk);
                                    den = den + (*denominator)(ii, jj, kk);
                                    nAvg = nAvg + 1;
                                }
                            }
                        }
                    }

                    (*numeratorAvg)(i, j, k) = num / nAvg;
                    (*denominatorAvg)(i, j, k) = den / nAvg;
                }
            }
        }
    }
    else
    {
        //! average in I direction
        if (averageDirection[0] != 0)
        {
            int nAvg = ied - ist + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int i = ist; i <= ied; ++ i)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*numeratorAvg)(i, j, k) = num / static_cast<RDouble>(nAvg);
                        (*denominatorAvg)(i, j, k) = den / static_cast<RDouble>(nAvg);
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }

        //! average in J direction
        if (averageDirection[1] != 0)
        {
            int nAvg = jed - jst + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int j = jst; j <= jed; ++ j)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*numeratorAvg)(i, j, k) = num / static_cast<RDouble>(nAvg);
                        (*denominatorAvg)(i, j, k) = den / static_cast<RDouble>(nAvg);
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }

        //! average in K direction
        if (averageDirection[2] != 0)
        {
            int nAvg = ked - kst + 1;

            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int k = kst; k <= ked; ++ k)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*numeratorAvg)(i, j, k) = num / static_cast<RDouble>(nAvg);
                        (*denominatorAvg)(i, j, k) = den / static_cast<RDouble>(nAvg);
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }
    }

    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    //isotropicConstant(I, J, K) = (*numeratorAvg)(I, J, K) / (*denominatorAvg)(I, J, K);
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble num = (*numeratorAvg)(i, j, k);
                RDouble den = (*denominatorAvg)(i, j, k);
                isotropicConstant(i, j, k) = num / den;
            }
        }
    }

    GhostCell3D(isotropicConstant, ni, nj, nk);

    delete numerator;
    delete denominator;
    delete numeratorAvg;
    delete denominatorAvg;
}

void LESSolverStruct::ComputePrandtlNumber(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *averageDirection = parameters->GetAverageDirection();

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));

    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));

    RDouble4D &leonardTemperatureFlux = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardTemperatureFlux"));
    RDouble4D &modeledTemperatureFlux = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledTemperatureFlux"));

    RDouble3D *numerator = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *denominator = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *numeratorAvg = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *denominatorAvg = new RDouble3D(I, J, K, fortranArray);

    ComputeLeonardTemperatureTerm(grid);

    ComputeModeledTemperatureTerm(grid);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble mTF0 = modeledTemperatureFlux(i, j, k, 0);
                RDouble mTF1 = modeledTemperatureFlux(i, j, k, 1);
                RDouble mTF2 = modeledTemperatureFlux(i, j, k, 2);

                RDouble lTF0 = leonardTemperatureFlux(i, j, k, 0);
                RDouble lTF1 = leonardTemperatureFlux(i, j, k, 1);
                RDouble lTF2 = leonardTemperatureFlux(i, j, k, 2);

                (*numerator)(i, j, k) = mTF0 * mTF0 + mTF1 * mTF1 + mTF2 * mTF2;
                (*denominator)(i, j, k) = lTF0 * mTF0 + lTF1 * mTF1 + lTF2 * mTF2;
            }
        }
    }

    (*numeratorAvg)(I, J, K) = (*numerator)(I, J, K);
    (*denominatorAvg)(I, J, K) = (*denominator)(I, J, K);

    if (averageDirection[0] == 0 && averageDirection[1] == 0 && averageDirection[2] == 0)
    {
        int averageWidth = parameters->GetAverageWidth();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    int nAvg = 0;

                    for (int ka = -averageWidth; ka <= averageWidth; ++ ka)
                    {
                        for (int ja = -averageWidth; ja <= averageWidth; ++ ja)
                        {
                            for (int ia = -averageWidth; ia <= averageWidth; ++ ia)
                            {
                                int kk = k + ka;
                                int jj = j + ja;
                                int ii = i + ia;

                                bool ifAverage = (ii >= ist && ii <= ied) && (jj >= jst && jj <= jed) && (kk >= kst && kk <= ked);

                                if (ifAverage)
                                {
                                    num = num + (*numerator)(ii, jj, kk);
                                    den = den + (*denominator)(ii, jj, kk);
                                    nAvg = nAvg + 1;
                                }
                            }
                        }
                    }

                    (*numeratorAvg)(i, j, k) = num / nAvg;
                    (*denominatorAvg)(i, j, k) = den / nAvg;
                }
            }
        }
    }
    else
    {
        //! average in I direction
        if (averageDirection[0] != 0)
        {
            int nAvg = ied - ist + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int i = ist; i <= ied; ++ i)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*numeratorAvg)(i, j, k) = num / nAvg;
                        (*denominatorAvg)(i, j, k) = den / nAvg;
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }

        //! average in J direction
        if (averageDirection[1] != 0)
        {
            int nAvg = jed - jst + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int j = jst; j <= jed; ++ j)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*numeratorAvg)(i, j, k) = num / nAvg;
                        (*denominatorAvg)(i, j, k) = den / nAvg;
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }

        //! average in K direction
        if (averageDirection[2] != 0)
        {
            int nAvg = ked - kst + 1;

            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int k = kst; k <= ked; ++ k)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*numeratorAvg)(i, j, k) = num / nAvg;
                        (*denominatorAvg)(i, j, k) = den / nAvg;
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }
    }

    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble cr = anisotropicConstant(i, j, k);

                RDouble KjKj = (*numeratorAvg)(i, j, k);
                RDouble NjNj = (*denominatorAvg)(i, j, k);

                RDouble prt = cr * KjKj / NjNj;
                turbulentPrandtlNumber(i, j, k) = prt;
            }
        }
    }

    GhostCell3D(turbulentPrandtlNumber, ni, nj, nk);

    //turbulentPrandtlNumber(I, J, K) = anisotropicConstant(I, J, K) * (*numeratorAvg)(I, J, K) / (*denominatorAvg)(I, J, K);

    delete numerator;
    delete denominator;
    delete numeratorAvg;
    delete denominatorAvg;
}

void LESSolverStruct::ComputeLeonardTemperatureTerm(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    RDouble4D &leonardTemperatureFlux = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardTemperatureFlux"));

    RDouble4D &testFilteredRhoUT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUT"));
    //RDouble3D &testFilteredDensity = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredDensity"));
    RDouble3D &testFilteredTemperature = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredTemperature"));
    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));

    ComputeTestFilteredRhoUT(grid);

    ComputeTestFilteredTemperature(grid);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble tFRhoUT = testFilteredRhoUT(i, j, k, 1);
                RDouble tFRhoVT = testFilteredRhoUT(i, j, k, 2);
                RDouble tFRhoWT = testFilteredRhoUT(i, j, k, 3);

                RDouble tFRho = testFilteredQ(i, j, k, 0);

                RDouble tFU = testFilteredQ(i, j, k, 1);
                RDouble tFV = testFilteredQ(i, j, k, 2);
                RDouble tFW = testFilteredQ(i, j, k, 3);

                RDouble tFT = testFilteredTemperature(i, j, k);

                leonardTemperatureFlux(i, j, k, 0) = tFRhoUT - tFRho * tFU * tFT;
                leonardTemperatureFlux(i, j, k, 1) = tFRhoVT - tFRho * tFV * tFT;
                leonardTemperatureFlux(i, j, k, 2) = tFRhoWT - tFRho * tFW * tFT;
            }
        }
    }
}

void LESSolverStruct::ComputeModeledTemperatureTerm(Grid *gridIn)
{
    using namespace IDX;
 
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Param_LESSolverStruct *parameters = GetControlParameters();
    RDouble testFilterScale = parameters->GetTestFilterScale();

    RDouble4D &modeledTemperatureFlux = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledTemperatureFlux"));

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));
    RDouble4D &testFilteredRhoStrainRateMagnitudeGradT = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoStrainRateMagnitudeGradT"));
    RDouble4D &testFilteredGradT = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredGradT"));
    RDouble4D &testFilteredStrainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredStrainRateTensor"));

    ComputeTestFilteredGradT(grid);

    ComputeTestFilteredRhoStrainRateMagnitudeGradT(grid);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble tFRSGradTx = testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 0);
                RDouble tFRSGradTy = testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 1);
                RDouble tFRSGradTz = testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 2);

                RDouble tFGradTx = testFilteredGradT(i, j, k, 0);
                RDouble tFGradTy = testFilteredGradT(i, j, k, 1);
                RDouble tFGradTz = testFilteredGradT(i, j, k, 2);

                RDouble tFRho = testFilteredQ(i, j, k, 0);
                RDouble tFS = testFilteredStrainRateTensor(i, j, k, 0);

                modeledTemperatureFlux(i, j, k, 0) = tFRSGradTx - testFilterScale * testFilterScale * tFRho * tFS * tFGradTx;
                modeledTemperatureFlux(i, j, k, 1) = tFRSGradTy - testFilterScale * testFilterScale * tFRho * tFS * tFGradTy;
                modeledTemperatureFlux(i, j, k, 2) = tFRSGradTz - testFilterScale * testFilterScale * tFRho * tFS * tFGradTz;
            }
        }
    }
}

void LESSolverStruct::ComputeTestFilteredRhoUT(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &temperatures = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("t"));

    RDouble4D &testFilteredRhoUT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUT"));

    testFilteredRhoUT(I, J, K, 1) = q(I, J, K, IR) * q(I, J, K, IU) * temperatures(I, J, K, ITT);
    testFilteredRhoUT(I, J, K, 2) = q(I, J, K, IR) * q(I, J, K, IV) * temperatures(I, J, K, ITT);
    testFilteredRhoUT(I, J, K, 3) = q(I, J, K, IR) * q(I, J, K, IW) * temperatures(I, J, K, ITT);

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    //! filter rhoUT in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 1; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredRhoUT(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredRhoUT(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredRhoUT(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter rhoUT in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 1; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredRhoUT(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredRhoUT(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredRhoUT(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter rhoUT in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 1; n <= 3; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredRhoUT(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredRhoUT(i, j, K, n);                
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredRhoUT(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredGradT(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &gradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradT"));
    RDouble4D &testFilteredGradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredGradT"));

    //RDouble4D &gradUVWTCellCenterX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    //RDouble4D &gradUVWTCellCenterY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    //RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    ComputeGradT(grid);

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    //! filter gradT in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 0; n <= 2; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredGradT(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredGradT(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredGradT(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    GhostCell3D(testFilteredGradT, ni, nj, nk, 3);

    //! filter gradT in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 1; n <= 2; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredGradT(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredGradT(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredGradT(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    GhostCell3D(testFilteredGradT, ni, nj, nk, 3);

    //! filter gradT in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 1; n <= 2; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredGradT(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredGradT(i, j, K, n);                
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredGradT(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    GhostCell3D(testFilteredGradT, ni, nj, nk, 3);

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredRhoStrainRateMagnitudeGradT(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &testFilteredRhoStrainRateMagnitudeGradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoStrainRateMagnitudeGradT"));
    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &strainRateTensor = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));
    RDouble4D &gradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradT"));

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble rho = q(i, j, k, IR);
                RDouble strainRateMagnitude = strainRateTensor(i, j, k, 0);
                RDouble dtdx = gradT(i, j, k, 0);
                RDouble dtdy = gradT(i, j, k, 1);
                RDouble dtdz = gradT(i, j, k, 2);

                testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 0) = rho * strainRateMagnitude * dtdx;
                testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 1) = rho * strainRateMagnitude * dtdy;
                testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 2) = rho * strainRateMagnitude * dtdz;
            }
        }
    }

    GhostCell3D(testFilteredRhoStrainRateMagnitudeGradT, ni, nj, nk, 3);

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    //! filter strainRateMagnitudeGradT in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 0; n <= 2; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredRhoStrainRateMagnitudeGradT(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredRhoStrainRateMagnitudeGradT(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredRhoStrainRateMagnitudeGradT(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    GhostCell3D(testFilteredRhoStrainRateMagnitudeGradT, ni, nj, nk, 3);

    //! filter strainRateMagnitudeGradT in J direction
    if (filterDirection[1] == 1) 
    {
        for (int n = 0; n <= 2; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredRhoStrainRateMagnitudeGradT(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredRhoStrainRateMagnitudeGradT(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredRhoStrainRateMagnitudeGradT(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    GhostCell3D(testFilteredRhoStrainRateMagnitudeGradT, ni, nj, nk, 3);

    //! filter strainRateMagnitudeGradT in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 0; n <= 2; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredRhoStrainRateMagnitudeGradT(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredRhoStrainRateMagnitudeGradT(i, j, K, n);                
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredRhoStrainRateMagnitudeGradT(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    GhostCell3D(testFilteredRhoStrainRateMagnitudeGradT, ni, nj, nk, 3);

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

/*
void LESSolverStruct::WALEOld(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble5D &subgridStress = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridStress"));
    //RDouble4D &subgridEnergy = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("kSGS"));
    RDouble5D &subgridHeatFlux = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridHeatFlux"));

    RDouble5D &gradPrimtiveVarFaceXofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceXofLES"));
    RDouble5D &gradPrimtiveVarFaceYofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceYofLES"));
    RDouble5D &gradPrimtiveVarFaceZofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceZofLES"));

    RDouble5D &gradtx = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceXofLES"));
    RDouble5D &gradty = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceYofLES"));
    RDouble5D &gradtz = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceZofLES"));

    RDouble4D &faceLengthScale = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("faceLengthScale"));
    RDouble4D &wallFunction = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("wallFunction"));

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble refReNumber = parameters->GetRefReNumber();

    RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

    int turbViscousCutType = parameters->GetTurbViscousCutType();

    RDouble waleConstant = parameters->GetWaleConstant();

    RDouble isotropicConstant = parameters->GetIsotropicConstant();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    RDouble refGama = parameters->GetRefGama();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);

    double vistmax  = 0.0;
    double vistmin  = 1.0e30;
    double logsmall = log(SMALL);

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    int imin = 1;
    int jmin = 1;
    int kmin = 1;

    //ComputeSubgridScaleEnergy(grid);

    int nDim = GetDim();

    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int ist, ied, jst, jed, kst, ked;
        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble deltaBar = faceLengthScale(il ,jl ,kl, iSurface);

                    RDouble dudx = gradPrimtiveVarFaceXofLES(il, jl, kl, IU, iSurface);
                    RDouble dudy = gradPrimtiveVarFaceYofLES(il, jl, kl, IU, iSurface);
                    RDouble dudz = gradPrimtiveVarFaceZofLES(il, jl, kl, IU, iSurface);

                    RDouble dvdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IV, iSurface);

                    RDouble dwdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IW, iSurface);

                    RDouble dtdx = gradtx(il, jl, kl, ITT, iSurface);
                    RDouble dtdy = gradty(il, jl, kl, ITT, iSurface);
                    RDouble dtdz = gradtz(il, jl, kl, ITT, iSurface);

                    RDouble s11 = dudx;
                    RDouble s22 = dvdy;
                    RDouble s33 = dwdz;
                    RDouble s12 = 0.5 * (dudy + dvdx);
                    RDouble s13 = 0.5 * (dudz + dwdx);
                    RDouble s23 = 0.5 * (dvdz + dwdy);
                    RDouble omega12 = 0.5 * (dudy - dvdx);
                    RDouble omega13 = 0.5 * (dudz - dwdx);
                    RDouble omega23 = 0.5 * (dvdz - dwdy);

                    RDouble skk = s11 + s22 + s33;
                    RDouble sij2 = s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
                    RDouble omegaij2 = 2.0 * (omega12 * omega12 + omega13 * omega13 + omega23 * omega23);

                    RDouble ss = 1.0/3.0 * (sij2 - omegaij2);
                    RDouble sd11 = s11 * s11 + s12 * s12 + s13 * s13 - omega12 * omega12 - omega13 * omega13 - ss;
                    RDouble sd22 = s12 * s12 + s22 * s22 + s23 * s23 - omega12 * omega12 - omega23 * omega23 - ss;
                    RDouble sd33 = s13 * s13 + s23 * s23 + s33 * s33 - omega13 * omega13 - omega23 * omega23 - ss;
                    RDouble sd12 = s11 * s12 + s12 * s22 + s13 * s23 - omega13 * omega23;
                    RDouble sd13 = s11 * s13 + s12 * s23 + s13 * s33 + omega12 * omega23;
                    RDouble sd23 = s12 * s13 + s22 * s23 + s23 * s33 - omega12 * omega13;

                    RDouble sdij2 = sd11 * sd11 + sd22 * sd22 + sd33 * sd33 + 2.0 * (sd12 * sd12 + sd13 * sd13 + sd23 * sd23);

                    RDouble rho = 0.5 * (q(i, j, k, IR) + q(il, jl, kl, IR));
                    RDouble mul = 0.5 * (visl(i, j, k) + visl(il, jl, kl));
                    RDouble dampingFunction = wallFunction(il, jl, kl, iSurface);

                    ss = (waleConstant * deltaBar * dampingFunction) * (waleConstant * deltaBar * dampingFunction) * refReNumber;
            
                    RDouble mut = rho * ss * pow(sdij2, 1.5) / (pow(sij2, 2.5) + pow(sdij2, 1.25));

                    if (turbViscousCutType == 1)
                    {
                        mut = MAX(mut - mul , 0.0);
                    }
                    else if (turbViscousCutType == 2)
                    {
                        mut = MAX(mut, 0.0);
                    }
                    else if (turbViscousCutType == 3)
                    {
                        mut = MAX(mut, -mul);
                    }

                    RDouble kSGS = 2.0 * isotropicConstant * (mut / deltaBar / dampingFunction) * (mut / deltaBar / dampingFunction) / rho / refReNumber;
                 
                    subgridStress(il, jl, kl, 1, iSurface) = -2.0 * mut * (s11 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    subgridStress(il, jl, kl, 2, iSurface) = -2.0 * mut * (s22 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    subgridStress(il, jl, kl, 3, iSurface) = -2.0 * mut * (s33 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    subgridStress(il, jl, kl, 4, iSurface) = -2.0 * mut * s12;
                    subgridStress(il, jl, kl, 5, iSurface) = -2.0 * mut * s13;
                    subgridStress(il, jl, kl, 6, iSurface) = -2.0 * mut * s23;

                    subgridHeatFlux(il, jl, kl, 1, iSurface) = -mut * specificHeatAtConstantPressure * oPrandtlTurbulence * dtdx;
                    subgridHeatFlux(il, jl, kl, 2, iSurface) = -mut * specificHeatAtConstantPressure * oPrandtlTurbulence * dtdy;
                    subgridHeatFlux(il, jl, kl, 3, iSurface) = -mut * specificHeatAtConstantPressure * oPrandtlTurbulence * dtdz;

                    if (vistmax < mut)
                    {
                        vistmax = mut;
                        imax = i;
                        jmax = j;
                        kmax = k;
                    }

                    if (vistmin > mut)
                    {
                        vistmin = mut;
                        imin = i;
                        jmin = j;
                        kmin = k;
                    }
                }
            }
        }
    }

    double vist_max = vistmax;
    double vist_min = vistmin;

    grid->UpdateData("vist_max", &vist_max, PHDOUBLE, 1);
    grid->UpdateData("vist_min", &vist_min, PHDOUBLE, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorVistmax)
        {
            cout << "vistmax = " << vistmax << " " << imax << " " << jmax << " " << kmax << "\n";
            cout << "vistmin = " << vistmin << " " << imin << " " << jmin << " " << kmin << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif

}
*/

void LESSolverStruct::WALE(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D &cellLengthScale = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cellLengthScale"));
    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    Param_LESSolverStruct *parameters = GetControlParameters();

    int turbViscousCutType = parameters->GetTurbViscousCutType();

    RDouble waleConstant = parameters->GetWaleConstant();

    RDouble iConstant = parameters->GetIsotropicConstant();

    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));
    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    double vistmax  = 0.0;
    double vistmin  = 1.0e30;
    double logsmall = log(SMALL);

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    int imin = 1;
    int jmin = 1;
    int kmin = 1;

    //ComputeSubgridScaleEnergy(grid);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {                
                RDouble deltaBar = cellLengthScale(i, j, k);

                RDouble dampingFunction = wallFunction(i, j, k);

                RDouble cw = (waleConstant * deltaBar * dampingFunction) * (waleConstant * deltaBar * dampingFunction);

                RDouble ci = iConstant * (deltaBar * dampingFunction) * (deltaBar * dampingFunction);

                anisotropicConstant(i, j, k) = cw;

                isotropicConstant(i, j, k) = ci;
            }
        }
    }

    RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();

    turbulentPrandtlNumber = prandtlTurbulence;
}

/*
void LESSolverStruct::ComputeSubgridScaleEnergy(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range IFACE, JFACE, KFACE;
    GetRange(ni, nj, nk, 0, 0, IFACE, JFACE, KFACE);

    Range Direction(1, 3);

    RDouble4D &faceLengthScale = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("faceLengthScale"));
    RDouble4D &wallFunction = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble4D &subgridScaleEnergy = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("kSGS"));
    RDouble3D &isotropicConstant = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));

    RDouble5D &gradPrimtiveVarFaceXofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceXofLES"));
    RDouble5D &gradPrimtiveVarFaceYofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceYofLES"));
    RDouble5D &gradPrimtiveVarFaceZofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceZofLES"));

    Param_LESSolverStruct * parameters = GetControlParameters();

    RDouble refReNumber  = parameters->GetRefReNumber();

    int subgridScaleModel = parameters->GetSubgridScaleModel();
    if (subgridScaleModel == 1 || subgridScaleModel == 3)
    {
        isotropicConstant(I, J, K) = parameters->GetIsotropicConstant();
    }
    else if (subgridScaleModel == 2)
    {
        ComputeIsotropicConstant(grid);
        faceLengthScale(IFACE, JFACE, KFACE, Direction) = 1.0;
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }

    int nDim = GetDim();

    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int ist, ied, jst, jed, kst, ked;
        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble deltaBar = faceLengthScale(il ,jl ,kl, iSurface);

                    RDouble dudx = gradPrimtiveVarFaceXofLES(il, jl, kl, IU, iSurface);
                    RDouble dudy = gradPrimtiveVarFaceYofLES(il, jl, kl, IU, iSurface);
                    RDouble dudz = gradPrimtiveVarFaceZofLES(il, jl, kl, IU, iSurface);

                    RDouble dvdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IV, iSurface);

                    RDouble dwdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IW, iSurface);

                    RDouble s11 = dudx;
                    RDouble s22 = dvdy;
                    RDouble s33 = dwdz;
                    RDouble s12 = half * (dudy + dvdx);
                    RDouble s13 = half * (dudz + dwdx);
                    RDouble s23 = half * (dvdz + dwdy);

                    RDouble sijsij = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                    RDouble ss = sqrt(two * sijsij);

                    RDouble rho = 0.5 * (q(i, j, k, IR) + q(il, jl, kl, IR));

                    RDouble dampingFunction = wallFunction(il, jl, kl, iSurface);

                    RDouble ci = 0.5 * (isotropicConstant(i, j, k) + isotropicConstant(il, jl, kl));

                    subgridScaleEnergy(il, jl, kl, iSurface) = 2.0 * ci * (deltaBar * dampingFunction) * (deltaBar * dampingFunction) * rho * ss * ss;
                }
            }
        }

    }
}
*/

void LESSolverStruct::filter( const int & nn
                            , RFloat1D * fu1, RFloat1D * fu2, RFloat1D * fu3
                            , RFloat1D * fuifuj11, RFloat1D * fuifuj22, RFloat1D * fuifuj33, RFloat1D * fuifuj12, RFloat1D * fuifuj13, RFloat1D * fuifuj23
                            , RFloat1D * fs11    , RFloat1D * fs22    , RFloat1D * fs33    , RFloat1D * fs12    , RFloat1D * fs13    , RFloat1D * fs23    
                            , RFloat1D * fsfs11  , RFloat1D * fsfs22  , RFloat1D * fsfs33  , RFloat1D * fsfs12  , RFloat1D * fsfs13  , RFloat1D * fsfs23 )
{
    //--- get TF(ui)
    TestFilter(nn, fu1);
    TestFilter(nn, fu2);
    TestFilter(nn, fu3);
    //--- get T( F(ui)F(uj) )
    TestFilter(nn, fuifuj11);
    TestFilter(nn, fuifuj22);
    TestFilter(nn, fuifuj33);
    TestFilter(nn, fuifuj12);
    TestFilter(nn, fuifuj13);
    TestFilter(nn, fuifuj23);
    //--- get T( F(Sij) )
    TestFilter(nn, fs11);
    TestFilter(nn, fs22);
    TestFilter(nn, fs33);
    TestFilter(nn, fs12);
    TestFilter(nn, fs13);
    TestFilter(nn, fs23);
    //--- get T( F(|S|)F(Sij) )
    TestFilter(nn, fsfs11);
    TestFilter(nn, fsfs22);
    TestFilter(nn, fsfs33);
    TestFilter(nn, fsfs12);
    TestFilter(nn, fsfs13);
    TestFilter(nn, fsfs23);
}

void LESSolverStruct::TestFilter(const int &nn, RFloat1D *fu)
{
    RFloat1D *tfu;	
    Range II(1,nn);
    tfu = new RFloat1D(II, fortranArray);
    (* tfu) = 0.0;

    Param_LESSolverStruct * parameters = GetControlParameters();

    int testFilterWidth = parameters->GetTestFilterWidth();

    RDouble * testFilterCoef = parameters->GetTestFilterCoef();

    for (int i = 1; i <= (nn-1); ++ i)
    {
        (*tfu)(i) = 0.0;
        for (int nw = -testFilterWidth; nw <= testFilterWidth; ++ nw)
        {
            RDouble coef = testFilterCoef[nw + testFilterWidth];
            RDouble f = (*fu)(i+nw);

            (*tfu)(i) = (*tfu)(i) + coef * f;
        }
    }

    for (int i = 1; i <= (nn-1); ++ i)
    {
        (*fu)(i) = (*tfu)(i);
    }

    delete tfu;
}

void LESSolverStruct::TestFilter(const int &n0, const int &nn, RDouble1D *&tf, RDouble1D *&f)
{
    Param_LESSolverStruct * parameters = GetControlParameters();

    int testFilterWidth = parameters->GetTestFilterWidth();

    RDouble * testFilterCoef = parameters->GetTestFilterCoef();

    for (int i = n0; i <= nn; ++ i)
    {
        RDouble sum = 0.0;
        for (int nw = -testFilterWidth; nw <= testFilterWidth; ++ nw)
        {
            RDouble coef = testFilterCoef[nw + testFilterWidth];
            RDouble x = (*f)(i + nw);
            sum = sum + coef * x;

            (*tf)(i) = sum;
        }
    }
}

void LESSolverStruct::TestFilter(const int &n0, const int &nn, RDouble *tf, RDouble *f)
{
    Param_LESSolverStruct * parameters = GetControlParameters();

    int testFilterWidth = parameters->GetTestFilterWidth();

    RDouble * testFilterCoef = parameters->GetTestFilterCoef();

    for (int i = n0; i <= nn; ++ i)
    {
        tf[i] = 0.0;
        for (int nw = -testFilterWidth; nw <= testFilterWidth; ++ nw)
        {
            tf[i] = tf[i] + testFilterCoef[nw + testFilterWidth] * f[i+nw];
        }
    }
}

void LESSolverStruct::Boundary(Grid *grid_in)
{
}

void LESSolverStruct::GetDeltaBar(const int deltaFunctionType, const RDouble deltai, const RDouble deltaj, const RDouble deltak, RDouble &deltaBar)
{
    if (deltaFunctionType == 1)
    {
        deltaBar = MAX(MAX(deltai, deltaj), deltak);
    }
    else if (deltaFunctionType == 2)
    {
        deltaBar = pow((deltaj * deltak * deltai), third);
    }
    else if (deltaFunctionType == 3)
    {
        GetFilterWidthOfScotti(deltai, deltaj, deltak, deltaBar);
    }
    else
    {
        TK_Exit::UnexpectedVarValue( "deltaFunctionType", deltaFunctionType);
    }
}

void LESSolverStruct::GetFilterWidthOfScotti(const RDouble deltai, const RDouble deltaj, const RDouble deltak, RDouble &deltaBar)
{
    RDouble deltaMax = MAX(MAX(deltai, deltaj), deltak);
    RDouble a1 = deltai / deltaMax;
    RDouble a2 = deltaj / deltaMax;
    RDouble a3 = deltak / deltaMax;

    if (a2 < a1)
    {
        RDouble atmp = a1;
        a1 = a2;
        a2 = atmp;
    }
    if (a3 < a1)
    {
        RDouble atmp = a1;
        a1 = a3;
        a3 = atmp;
    }
    if (a3 < a2)
    {
        RDouble atmp = a2;
        a2 = a3;
        a3 = atmp;
    }

    deltaBar = cosh(sqrt(4.0 / 27.0) * (log(a1) * log(a1) - log(a1) * log(a2) + log(a2) * log(a2)));
}

}