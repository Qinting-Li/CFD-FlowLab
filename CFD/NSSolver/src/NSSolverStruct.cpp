#include "NSSolverStruct.h"
#include "Glb_Dimension.h"
#include "Geo_StructBC.h"
#include "MultiGridOperation.h"
#include "GradientOperation.h"
#include "Force.h"
#include "OversetGridFactory.h"
#include "TK_Exit.h"
#include "Gradient.h"
#include "FieldProxy.h"
#include "IO_FileName.h"
#include "Geo_UnstructGrid.h"
#include "Param_NSSolverStruct.h"
#include "Math_Limiter.h"
#include "Post_ForceMoment.h"
#include "Residual.h"
#include "Limiter.h"
#include "Jacobian.h"
#include "TK_Log.h"

#ifdef USE_INHOUSE_VERSION
#include "BLUSGS.h"
#endif

#ifdef USE_HYFLOW_CODE
#include "RealGas.h"
#include "EquilibriumGas.h"
#else
#include "Gas.h"
#endif
#include "TK_Time.h"
using namespace std;

namespace PHSPACE
{
    using namespace GAS_SPACE;

    NSSolverStruct::NSSolverStruct()
    {
        for (int i = 0; i < 5; ++ i)
        {
            gradientCellCenterUVWT[i] = 0;
        }
        nWallBC = 1; //! The number of surface boundary conditions.
        nMaxSurfaceCell = 1; //! The number of surface cells.
        nIonizationFlag = 0; //! The flag of chemical reactions.
    }

    NSSolverStruct::~NSSolverStruct()
    {
        DeAllocateGlobalVariables();
        FreeControlParameters();
    }

    void NSSolverStruct::ReadParameter()
    {
        Param_NSSolverStruct *parameters = GetControlParameters();

        string structSchemeName = parameters->GetInviscidSchemeName();
        int structScheme = GetSchemeID(structSchemeName);
        GlobalDataBase::UpdateData("str_scheme", &structScheme, PHINT, 1);

        if (structScheme == ISCHEME_ROE)
        {
            parameters->SetEntropyFixCoefficients();

            Grid *grid = GetGrid();
            int zoneID = grid->GetZoneID();

            RDouble entropyAcoustic = parameters->GetRoeEntropyFixCoef1();

            //! The output of Entropy fix coefficient for ROE scheme.
            if (zoneID == 0)  //! When grid->GetZoneID() = 0 and GetCurrentProcessorID() = 0.
            {
                PrintToWindow("Entropy fix coefficient for ROE scheme: ", entropyAcoustic, "\n");
            }
        }

        string structLimiterName = GlobalDataBase::GetStrParaFromDB("str_limiter_name");
        int structLimiter = GetLimiterID(structLimiterName);
        GlobalDataBase::UpdateData("str_limiter", &structLimiter, PHINT, 1);

        string unstructLimiterName = GlobalDataBase::GetStrParaFromDB("uns_limiter_name");
        int unstructLimiter = GetLimiterID(unstructLimiterName);
        GlobalDataBase::UpdateData("uns_limiter", &unstructLimiter, PHINT, 1);
    }

    void NSSolverStruct::AllocateGlobalVar(Grid *gridIn)
    {
        CFDSolver::AllocateGlobalVar(gridIn);

        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Range I, J, K;
        GetRange(ni, nj, nk, -2, 1, I, J, K);

        Range IFACE, JFACE, KFACE;
        GetRange(ni, nj, nk, -2, 2, IFACE, JFACE, KFACE);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int isWennScheme = parameters->GetWennSchemeFlag();

        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
        int nSourceTerm = nEquation - nNSEquation;    //! The number of source term.
        int nDim = GetDim();
        const int nUVWTIndex = nTemperatureModel + 3;

        int nLayer = GetNumberOfGhostCellLayers();

        if (isWennScheme == 1)
        {
            GetRange(ni, nj, nk, -3, 2, I, J, K);
        }

        Range M(0, nEquation - 1);
        Range N(0, nUVWTIndex - 1);
        Range T(0, nTemperatureModel - 1);
        Range D(0, nDim - 1);

        RDouble4D *q = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *residual = new RDouble4D(I, J, K, M, fortranArray);
        RDouble3D *dt = new RDouble3D(I, J, K, fortranArray);
        RDouble3D *gama = new RDouble3D(I, J, K, fortranArray);
        RDouble4D *t = new RDouble4D(I, J, K, T, fortranArray);
        RDouble4D *visSpectralRadius = new RDouble4D(I, J, K, D, fortranArray);
        RDouble4D *invSpectralRadius = new RDouble4D(I, J, K, D, fortranArray);
        RDouble4D *diagonal = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *diagonal0 = new RDouble4D(I, J, K, M, fortranArray);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, nLayer);


        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    for (int it = 0; it < nTemperatureModel; ++ it)
                    {
                        (*t)(i, j, k, it) = 0.0;
                    }
                }
            }
        }

        grid->UpdateDataPtr("q", q);    //! Scalar flow field variable(rho/u/v/w/p).
        grid->UpdateDataPtr("res", residual);    //! Residual or right-hand side.
        grid->UpdateDataPtr("dt", dt);    //! Time step.
        grid->UpdateDataPtr("gama", gama);    //! Ratio of specific heat coefficients at constant pressure and volume.
        grid->UpdateDataPtr("t", t);    //! Static temperature.
        grid->UpdateDataPtr("visSpectralRadius", visSpectralRadius);    //! Viscous spectral radius.
        grid->UpdateDataPtr("invSpectralRadius", invSpectralRadius);    //! Inviscid spectral radius.
        grid->UpdateDataPtr("diagonal", diagonal);    //! Sacrificing space to improve efficiency.
        grid->UpdateDataPtr("diagonal0", diagonal0);    //! Sacrificing space to improve efficiency.

        //! Allocate memories.
        int timeIntegration = parameters->GetTimeIntegration();
        if (timeIntegration == 5)
        {
            RDouble5D *diagonalMatrice = new RDouble5D(I, J, K, M, M, fortranArray);
            grid->UpdateDataPtr("diagonalMatrice", diagonalMatrice);
        }
        if (nSourceTerm > 0)
        {
            //! To store the Jacobian matrix of source term.
            RDouble4D *sourceJacobian = new RDouble4D(I, J, K, Range(0, nSourceTerm * nSourceTerm - 1), fortranArray);
            RDouble4D *sourceChemical = new RDouble4D(I, J, K, Range(0, nSourceTerm - 1), fortranArray);
            RDouble5D *sourceDerivatives = new RDouble5D(I, J, K, Range(0, nSourceTerm - 1), Range(0, nSourceTerm - 1), fortranArray);
            grid->UpdateDataPtr("sourceJacobian", sourceJacobian);
            grid->UpdateDataPtr("sourceChemical", sourceChemical);
            grid->UpdateDataPtr("sourceDerivatives", sourceDerivatives);

            //! To store the partial derivatives of source terms in the main diagonal position.
            RDouble4D *srs = new RDouble4D(I, J, K, Range(0, nSourceTerm - 1), fortranArray);
            grid->UpdateDataPtr("srs", srs);
        }

        int isUnsteady = parameters->GetIsUnsteady();
        if (isUnsteady)
        {
            RDouble4D *q_unsteady_n1 = new RDouble4D(I, J, K, M, fortranArray);
            RDouble4D *q_unsteady_n2 = new RDouble4D(I, J, K, M, fortranArray);
            RDouble4D *res_unsteady_n1 = new RDouble4D(I, J, K, M, fortranArray);
            RDouble4D *res_unsteady_n2 = new RDouble4D(I, J, K, M, fortranArray);

            //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dual-time method.
            RDouble4D *res_unsteady_tmp = new RDouble4D(I, J, K, M, fortranArray);

            grid->UpdateDataPtr("q_unsteady_n1", q_unsteady_n1);
            grid->UpdateDataPtr("q_unsteady_n2", q_unsteady_n2);
            grid->UpdateDataPtr("res_unsteady_n1", res_unsteady_n1);
            grid->UpdateDataPtr("res_unsteady_n2", res_unsteady_n2);
            grid->UpdateDataPtr("res_unsteady_tmp", res_unsteady_tmp);

            //! Statistical variables for unsteady simulation.
            RDouble4D *qAverage = new RDouble4D(I, J, K, M, fortranArray);
            *qAverage = 0.0;
            grid->UpdateDataPtr("qAverage", qAverage);

            //! Statistical Reynolds stress for unsteady simulation.
            RDouble4D *tauAverage = new RDouble4D(I, J, K, Range(0, 5), fortranArray);
            *tauAverage = 0.0;
            grid->UpdateDataPtr("tauAverage", tauAverage);

            RDouble4D *q2Average = new RDouble4D(I, J, K, Range(0, 5), fortranArray);
            *q2Average = 0.0;
            grid->UpdateDataPtr("q2Average", q2Average);
        }

        RDouble3D *rtem = new RDouble3D(I, J, K, fortranArray);
        *rtem = 0.0;
        grid->UpdateDataPtr("rtem", rtem);    //! Pressure factor.

        //RDouble3D &rtem_initial = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("rtem"));
        //rtem_initial = 0.0;
        //! rtem is used in timestep(), it should be initialed.

        int viscousType = parameters->GetViscousType();    //! Viscous model.
        if (viscousType > INVISCID)
        {
            RDouble3D *viscousLaminar = new RDouble3D(I, J, K, fortranArray);
            RDouble3D *viscousTurbulence = new RDouble3D(I, J, K, fortranArray);
            RDouble4D *heatTransferCoeff = new RDouble4D(I, J, K, T, fortranArray);
            grid->UpdateDataPtr("heatTransferCoeff", heatTransferCoeff);
            grid->UpdateDataPtr("visl", viscousLaminar);    //! Laminar viscous coefficient.
            grid->UpdateDataPtr("vist", viscousTurbulence);    //! Turbulence viscous coefficient.
            *viscousLaminar = 1.0;
            *viscousTurbulence = 0.0;
            *heatTransferCoeff = 0.0;
        }

        int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
        if (iLES == LES_SOLVER)
        {
            RDouble3D *subgridScaleEnergy = new RDouble3D(I, J, K, fortranArray);
            grid->UpdateDataPtr("subgridScaleEnergy", subgridScaleEnergy);
            *subgridScaleEnergy = 0.0;
        }
        if (nSourceTerm > 0 && viscousType > INVISCID) //Save the heat conductivity and mass diffusion coefficient.
        {
            RDouble4D *heatConductivity = new RDouble4D(I, J, K, T, fortranArray);
            RDouble4D *massDiffusionCoef = new RDouble4D(I, J, K, Range(0, nSourceTerm - 1), fortranArray);
            grid->UpdateDataPtr("lambda", heatConductivity);    //! Thermal conductivity coefficient.
            grid->UpdateDataPtr("rhoD", massDiffusionCoef);    //! Mass Diffusion coefficient.
        }

        RDouble4D *gradPrimtiveVarFaceX = new RDouble4D(IFACE, JFACE, KFACE, M, fortranArray);
        RDouble4D *gradPrimtiveVarFaceY = new RDouble4D(IFACE, JFACE, KFACE, M, fortranArray);
        RDouble4D *gradPrimtiveVarFaceZ = new RDouble4D(IFACE, JFACE, KFACE, M, fortranArray);
        grid->UpdateDataPtr("gradPrimtiveVarFaceX", gradPrimtiveVarFaceX);    //! Gradient of scalar flow field variable at face, for x direction.
        grid->UpdateDataPtr("gradPrimtiveVarFaceY", gradPrimtiveVarFaceY);    //! Gradient of scalar flow field variable at face, for y direction.
        grid->UpdateDataPtr("gradPrimtiveVarFaceZ", gradPrimtiveVarFaceZ);    //! Gradient of scalar flow field variable at face, for z direction.

        RDouble4D *gradUVWTCellCenterX = new RDouble4D(I, J, K, N, fortranArray);
        RDouble4D *gradUVWTCellCenterY = new RDouble4D(I, J, K, N, fortranArray);
        RDouble4D *gradUVWTCellCenterZ = new RDouble4D(I, J, K, N, fortranArray);
        grid->UpdateDataPtr("gradUVWTCellCenterX", gradUVWTCellCenterX);    //! Gradient of velocity components(u/v/w) and temperature at cell center in x-direction.
        grid->UpdateDataPtr("gradUVWTCellCenterY", gradUVWTCellCenterY);    //! Gradient of velocity components(u/v/w) and temperature at cell center in y-direction.
        grid->UpdateDataPtr("gradUVWTCellCenterZ", gradUVWTCellCenterZ);    //! Gradient of velocity components(u/v/w) and temperature at cell center in z-direction.

        //! Temporary variable for periodic boundary condition.
        RDouble4D *rotNSgradValueX = new RDouble4D(I, J, K, N, fortranArray);
        RDouble4D *rotNSgradValueY = new RDouble4D(I, J, K, N, fortranArray);
        RDouble4D *rotNSgradValueZ = new RDouble4D(I, J, K, N, fortranArray);
        grid->UpdateDataPtr("rotNSgradValueX", rotNSgradValueX);
        grid->UpdateDataPtr("rotNSgradValueY", rotNSgradValueY);
        grid->UpdateDataPtr("rotNSgradValueZ", rotNSgradValueZ);

        RDouble4D *gradTemperatureFaceX = new RDouble4D(IFACE, JFACE, KFACE, T, fortranArray);
        RDouble4D *gradTemperatureFaceY = new RDouble4D(IFACE, JFACE, KFACE, T, fortranArray);
        RDouble4D *gradTemperatureFaceZ = new RDouble4D(IFACE, JFACE, KFACE, T, fortranArray);
        grid->UpdateDataPtr("gradTemperatureFaceX", gradTemperatureFaceX);    //! Gradient of static temperature at face, for x direction.
        grid->UpdateDataPtr("gradTemperatureFaceY", gradTemperatureFaceY);    //! Gradient of static temperature at face, for y direction.
        grid->UpdateDataPtr("gradTemperatureFaceZ", gradTemperatureFaceZ);    //! Gradient of static temperature at face, for z direction.

        RDouble4D *ql1 = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *qr1 = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *ql2 = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *qr2 = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *ql3 = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *qr3 = new RDouble4D(I, J, K, M, fortranArray);
        grid->UpdateDataPtr("ql1", ql1);
        grid->UpdateDataPtr("qr1", qr1);
        grid->UpdateDataPtr("ql2", ql2);
        grid->UpdateDataPtr("qr2", qr2);
        grid->UpdateDataPtr("ql3", ql3);
        grid->UpdateDataPtr("qr3", qr3);

        int nSpeciesNumber = parameters->GetNumberOfSpecies();
        int nTotalNumber = 16;
        Int1D *speciesOrder = new Int1D(Range(0, nTotalNumber - 1), fortranArray);
        grid->UpdateDataPtr("speciesOrder", speciesOrder);
        for (int n = 0; n < nTotalNumber; ++ n)
        {
            (*speciesOrder)(n) = -1;
        }
        int iCodeBranch = GlobalDataBase::GetIntParaFromDB("iCodeBranch");
        if (iCodeBranch > 0)
        {
            RDouble *maxResidualVariation = new RDouble[nEquation];
            for (int m = 0; m < nEquation; ++m)
            {
                maxResidualVariation[m] = 0.0;
            }
            grid->UpdateDataPtr("maxResidualVariation", maxResidualVariation);
        }

        if (nChemical > 0 && nSpeciesNumber > 0)
        {
            string *nameList = gas->GetNameOfSpecies();
            SpeciesNameToInteger(nSpeciesNumber, nameList, speciesOrder);

            RDouble4D *gradSpeciesCellCenterX = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);
            RDouble4D *gradSpeciesCellCenterY = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);
            RDouble4D *gradSpeciesCellCenterZ = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);
            grid->UpdateDataPtr("dcdx", gradSpeciesCellCenterX);    //! Gradient of chemical species at cell center in x-direction.
            grid->UpdateDataPtr("dcdy", gradSpeciesCellCenterY);    //! Gradient of chemical species at cell center in y-direction.
            grid->UpdateDataPtr("dcdz", gradSpeciesCellCenterZ);    //! Gradient of chemical species at cell center in z-direction.

            RDouble4D *speciesEnthalpy = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);
            grid->UpdateDataPtr("hSpecies", speciesEnthalpy);    //! Enthalpy of chemical species.

            RDouble3D *totalEnergy = new RDouble3D(I, J, K, fortranArray);
            RDouble3D *totalEnthalpy = new RDouble3D(I, J, K, fortranArray);
            RDouble3D *totalCp = new RDouble3D(I, J, K, fortranArray);
            RDouble3D *totalCv = new RDouble3D(I, J, K, fortranArray);

            grid->UpdateDataPtr("totalEnergy", totalEnergy);    //! total energy.
            grid->UpdateDataPtr("totalEnthalpy", totalEnthalpy);    //! total Enthalpy.
            grid->UpdateDataPtr("totalCp", totalCp);    //! total Cp.
            grid->UpdateDataPtr("totalCv", totalCv);    //! total Cv.

            int nEnergyRecycle = parameters->GetnEnergyRecycle();
            if (nTemperatureModel > 1 && nEnergyRecycle == 2)
            {
                RDouble4D *temporaryEve = new RDouble4D(I, J, K, Range(1, MAX(1, nTemperatureModel - 1)), fortranArray);
                grid->UpdateDataPtr("temporaryEve", temporaryEve);    //! temporaryEve.
            }

            RDouble3D *totalCvtr = new RDouble3D(I, J, K, fortranArray);
            RDouble3D *totalCvv = new RDouble3D(I, J, K, fortranArray);
            RDouble3D *totalCve = new RDouble3D(I, J, K, fortranArray);
            RDouble4D *speciesCvs = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);
            RDouble4D *speciesCps = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);
            RDouble4D *speciesCvvs = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);
            RDouble4D *speciesCves = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);
            RDouble4D *speciesEtrs = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);
            RDouble4D *speciesEvs = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);
            RDouble4D *speciesEes = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);
            RDouble4D *speciesEs = new RDouble4D(I, J, K, Range(0, nSpeciesNumber - 1), fortranArray);


            grid->UpdateDataPtr("totalCvtr", totalCvtr);   //! total Cvtr.
            grid->UpdateDataPtr("totalCvv", totalCvv);    //! total Cvv.
            grid->UpdateDataPtr("totalCve", totalCve);    //! total Cve.
            grid->UpdateDataPtr("speciesCvs", speciesCvs);    //! Cvs for species.
            grid->UpdateDataPtr("speciesCps", speciesCps);    //! Cps for species.
            grid->UpdateDataPtr("speciesCvvs", speciesCvvs);    //! Cvvs for species.
            grid->UpdateDataPtr("speciesCves", speciesCves);    //! Cves for species.
            grid->UpdateDataPtr("speciesEtrs", speciesEtrs);    //! Etrs for species.
            grid->UpdateDataPtr("speciesEvs", speciesEvs);    //! Evs for species.
            grid->UpdateDataPtr("speciesEes", speciesEes);    //! Ees for species.
            grid->UpdateDataPtr("speciesEs", speciesEs);    //! Es for species.

            RDouble3D *tflag = new RDouble3D(I, J, K, fortranArray);
            *tflag = 0.0;

            grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd);
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        (*tflag)(i, j, k) = 1.0;
                    }
                }
            }
            grid->UpdateDataPtr("tflag", tflag);

#ifdef USE_HYFLOW_CODE
            modifyFlagForWall(grid);
#endif 

            grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        }

        //! Save the temperature on wall.
        GetSurfaceCellNumber(gridIn);
        RDouble2D *firstLayerHeight = NULL;
        if (nWallBC > 0)    //! The current grid includes wall boundary conditions.
        {
            firstLayerHeight = new RDouble2D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), fortranArray);
        }
        grid->UpdateDataPtr("firstLayerHeight", firstLayerHeight);
        //Compute the heights of cells on the first grid layer.
        ComputeFirstLayerGridHeight(gridIn);

        //Allocate the memories for temperatures on wall.
        if (nWallBC > 0)
        {
            RDouble2D *temperatureWall = new RDouble2D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), fortranArray);
            grid->UpdateDataPtr("surfaceTemperature", temperatureWall);
        }

        //Allocate the memories for catalytic wall.
        if (nChemical > 0 && nSpeciesNumber > 0)
        {
            RDouble3D *surfaceMassFraction = NULL;
            if (nWallBC > 0)
            {
                surfaceMassFraction = new RDouble3D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), Range(0, nSpeciesNumber - 1), fortranArray);
            }
            grid->UpdateDataPtr("surfaceMassFraction", surfaceMassFraction);
        }

        int isAblationWall = GlobalDataBase::GetIntParaFromDB("nAblation");
        int Injection = GlobalDataBase::GetIntParaFromDB("isInjection");
        if (isAblationWall > 0 && Injection > 0)
        {
            RDouble2D *injectionVelocity = NULL;
            if (nWallBC > 0)
            {
                injectionVelocity = new RDouble2D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), fortranArray);
            }

            for (int iWall = 0; iWall < nWallBC; ++ iWall)
            {
                for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
                {
                    (*injectionVelocity)(iWall, iCell) = 0.0;
                }
            }
            grid->UpdateDataPtr("injectionVelocity", injectionVelocity);
        }

        //! added by myk 202010, for conservation of flux
        systemgridtype = GlobalDataBase::GetIntParaFromDB("sys_gridtype");
        if (MIXGRID == systemgridtype)
        {
            RDouble4D *dqdx_cc_ruvwpt = new RDouble4D(I, J, K, Range(0, 5), fortranArray);
            RDouble4D *dqdy_cc_ruvwpt = new RDouble4D(I, J, K, Range(0, 5), fortranArray);
            RDouble4D *dqdz_cc_ruvwpt = new RDouble4D(I, J, K, Range(0, 5), fortranArray);
            grid->UpdateDataPtr("dqdx_cc_ruvwpt", dqdx_cc_ruvwpt);    //! Gradient of flow field variable(rho/u/v/w/p/t) at cell center in x-direction.
            grid->UpdateDataPtr("dqdy_cc_ruvwpt", dqdy_cc_ruvwpt);    //! Gradient of flow field variable(rho/u/v/w/p/t) at cell center in y-direction.
            grid->UpdateDataPtr("dqdz_cc_ruvwpt", dqdz_cc_ruvwpt);    //! Gradient of flow field variable(rho/u/v/w/p/t) at cell center in z-direction.

            RDouble4D *faceInviscidfluxI = new RDouble4D(I, J, K, M, fortranArray);
            RDouble4D *faceInviscidfluxJ = new RDouble4D(I, J, K, M, fortranArray);
            RDouble4D *faceInviscidfluxK = new RDouble4D(I, J, K, M, fortranArray);
            grid->UpdateDataPtr("faceInviscidfluxI", faceInviscidfluxI);    //! Inviscid flux at face, for I direction.
            grid->UpdateDataPtr("faceInviscidfluxJ", faceInviscidfluxJ);    //! Inviscid flux at face, for J direction.
            grid->UpdateDataPtr("faceInviscidfluxK", faceInviscidfluxK);    //! Inviscid flux at face, for K direction.

            if (viscousType > INVISCID)
            {
                RDouble4D *faceViscousfluxI = new RDouble4D(I, J, K, M, fortranArray);
                RDouble4D *faceViscousfluxJ = new RDouble4D(I, J, K, M, fortranArray);
                RDouble4D *faceViscousfluxK = new RDouble4D(I, J, K, M, fortranArray);
                grid->UpdateDataPtr("faceViscousfluxI", faceViscousfluxI);    //! Viscous flux at face, for I direction.
                grid->UpdateDataPtr("faceViscousfluxJ", faceViscousfluxJ);    //! Viscous flux at face, for J direction.
                grid->UpdateDataPtr("faceViscousfluxK", faceViscousfluxK);    //! Viscous flux at face, for K direction.
            }
        }

        //Allocate the memories for the slip parameters.
        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        if (nSlipBCModel > 0)
        {
            //The array includes variables such as slip temperatures(Tts, Tvs, Tes), slip velocity(Vx, Vy, Vz) and slip mass fractions of species.
            int nSlipVar = nSpeciesNumber + 6;
            RDouble3D *surfaceSlipVariables = NULL;
            if (nWallBC > 0)
            {
                surfaceSlipVariables = new RDouble3D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), Range(0, nSlipVar - 1), fortranArray);
            }
            grid->UpdateDataPtr("surfaceSlipVariables", surfaceSlipVariables);
        }

        int nRapidFlowfield = parameters->GetRapidFlowfieldMethod();
        if (nRapidFlowfield > 0)
        {
            RDouble3D *surfacePressure = NULL;
            if (nWallBC > 0)
            {
                surfacePressure = new RDouble3D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), Range(0, 1), fortranArray);
            }
            grid->UpdateDataPtr("surfacePressure", surfacePressure);
#ifdef USE_HYFLOW_CODE
            ComputeSurfacePressure(grid);
#endif
        }

        //To exam the surface heating change.
        int nSurfHeatMonitor = parameters->GetSurfaceHeatingMonitor();
        if (nSurfHeatMonitor > 0)
        {
            RDouble3D *surfaceHeatFlux = NULL;
            if (nWallBC > 0)
            {
                surfaceHeatFlux = new RDouble3D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), Range(0, 1), fortranArray);
                *surfaceHeatFlux = 0.0;
            }
            grid->UpdateDataPtr("surfaceHeatFlux", surfaceHeatFlux);
        }

        //initialization of parameters on surface.
        InitializeSurfaceParameters(gridIn);

        //Set the temperature on wall by reading data from file.
        InitializeSurfaceParametersFromFile(gridIn);
    }

    void NSSolverStruct::DeAllocateGlobalVar(Grid *gridIn)
    {
        CFDSolver::DeAllocateGlobalVar(gridIn);
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D *q = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D *residual = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble3D *dt = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));
        RDouble3D *gama = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));
        RDouble4D *t = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble4D *visSpectralRadius = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
        RDouble4D *invSpectralRadius = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
        RDouble4D *diagonal = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
        RDouble4D *diagonal0 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal0"));
        delete q; q = 0;
        delete residual; residual = 0;
        delete dt; dt = 0;
        delete gama; gama = 0;
        delete t; t = 0;
        delete visSpectralRadius; visSpectralRadius = 0;
        delete invSpectralRadius; invSpectralRadius = 0;
        delete diagonal; diagonal = 0;
        //Remove data pointer.
        grid->DeleteDataPtr("q");
        grid->DeleteDataPtr("res");
        grid->DeleteDataPtr("dt");
        grid->DeleteDataPtr("gama");
        grid->DeleteDataPtr("t");
        grid->DeleteDataPtr("visSpectralRadius");
        grid->DeleteDataPtr("invSpectralRadius");
        grid->DeleteDataPtr("diagonal");
        grid->DeleteDataPtr("diagonal0");

        Param_NSSolverStruct *parameters = GetControlParameters();
        int timeIntegration = parameters->GetTimeIntegration();

        if (timeIntegration == 5)
        {
            RDouble5D *diagonalMatrice = reinterpret_cast <RDouble5D *> (grid->GetDataPtr("diagonalMatrice"));
            delete diagonalMatrice; diagonalMatrice = 0;
            grid->DeleteDataPtr("diagonalMatrice");
        }

        int nEquation = GetNumberOfEquations();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nSourceTerm = nEquation - nNSEquation;     //! Number of source term.
        if (nSourceTerm > 0)
        {
            RDouble4D *sourceJacobian = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("sourceJacobian"));
            RDouble4D *sourceChemical = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("sourceChemical"));
            RDouble5D *sourceDerivatives = reinterpret_cast <RDouble5D *> (grid->GetDataPtr("sourceDerivatives"));
            delete sourceJacobian; sourceJacobian = 0;
            delete sourceChemical; sourceChemical = 0;
            delete sourceDerivatives; sourceDerivatives = 0;
            grid->DeleteDataPtr("sourceJacobian");
            grid->DeleteDataPtr("sourceChemical");
            grid->DeleteDataPtr("sourceDerivatives");

            RDouble4D *srs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("srs"));
            delete srs; srs = 0;
            grid->DeleteDataPtr("srs");
        }

        int isUnsteady = parameters->GetIsUnsteady();
        if (isUnsteady)
        {
            RDouble4D *q_unsteady_n1 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
            RDouble4D *q_unsteady_n2 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));
            RDouble4D *res_unsteady_n1 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n1"));
            RDouble4D *res_unsteady_n2 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n2"));

            //! Store the historical sum of InviscidFlux, ViscousFlux, et al., /n
            //! exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dual-time method.
            RDouble4D *res_unsteady_tmp = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_tmp"));

            delete q_unsteady_n1; q_unsteady_n1 = 0;
            delete q_unsteady_n2; q_unsteady_n2 = 0;
            delete res_unsteady_n1; res_unsteady_n1 = 0;
            delete res_unsteady_n2; res_unsteady_n2 = 0;
            delete res_unsteady_tmp; res_unsteady_tmp = 0;
            grid->DeleteDataPtr("q_unsteady_n1");
            grid->DeleteDataPtr("q_unsteady_n2");
            grid->DeleteDataPtr("res_unsteady_n1");
            grid->DeleteDataPtr("res_unsteady_n2");
            grid->DeleteDataPtr("res_unsteady_tmp");

            //! Statistical variables for unsteady simulation.
            int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
            if (ifStaticsFlowField > 0)
            {
                RDouble4D *qAverage = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage"));
                delete qAverage; qAverage = 0;
                grid->DeleteDataPtr("qAverage");
            }

            //! Statistical Reynolds stress for unsteady simulation.
            RDouble4D *tauAverage = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("tauAverage"));
            delete tauAverage; tauAverage = 0;
            grid->DeleteDataPtr("tauAverage");

            RDouble4D *q2Average = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q2Average"));
            delete q2Average; q2Average = 0;
            grid->DeleteDataPtr("q2Average");
        }

        RDouble3D *rtem = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("rtem"));
        delete rtem; rtem = 0;
        grid->DeleteDataPtr("rtem");

        int viscousType = parameters->GetViscousType();
        if (viscousType > INVISCID)
        {
            RDouble3D *viscousLaminar = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
            RDouble3D *viscousTurbulence = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
            RDouble4D *heatTransferCoeff = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("heatTransferCoeff"));
            delete viscousLaminar; viscousLaminar = 0;
            delete viscousTurbulence; viscousTurbulence = 0;
            delete heatTransferCoeff; heatTransferCoeff = 0;
            grid->DeleteDataPtr("visl");
            grid->DeleteDataPtr("vist");
            grid->DeleteDataPtr("heatTransferCoeff");
        }

        int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
        if (iLES == LES_SOLVER)
        {
            RDouble3D *subgridScaleEnergy = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));
            delete subgridScaleEnergy; subgridScaleEnergy = 0;
            grid->DeleteDataPtr("subgridScaleEnergy");
        }

        if (nSourceTerm > 0 && viscousType > INVISCID) //Save the heat conductivity and mass diffusion coefficient.
        {
            RDouble4D *heatConductivity = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("lambda"));
            RDouble4D *massDiffusionCoef = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rhoD"));
            delete heatConductivity; heatConductivity = 0;
            delete massDiffusionCoef; massDiffusionCoef = 0;
            grid->DeleteDataPtr("lambda");
            grid->DeleteDataPtr("rhoD");
        }

        RDouble4D *gradPrimtiveVarFaceX = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
        RDouble4D *gradPrimtiveVarFaceY = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
        RDouble4D *gradPrimtiveVarFaceZ = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));
        delete gradPrimtiveVarFaceX; gradPrimtiveVarFaceX = 0;
        delete gradPrimtiveVarFaceY; gradPrimtiveVarFaceY = 0;
        delete gradPrimtiveVarFaceZ; gradPrimtiveVarFaceZ = 0;
        grid->DeleteDataPtr("gradPrimtiveVarFaceX");
        grid->DeleteDataPtr("gradPrimtiveVarFaceY");
        grid->DeleteDataPtr("gradPrimtiveVarFaceZ");

        RDouble4D *gradTemperatureFaceX = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceX"));
        RDouble4D *gradTemperatureFaceY = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceY"));
        RDouble4D *gradTemperatureFaceZ = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceZ"));
        delete gradTemperatureFaceX; gradTemperatureFaceX = 0;
        delete gradTemperatureFaceY; gradTemperatureFaceY = 0;
        delete gradTemperatureFaceZ; gradTemperatureFaceZ = 0;
        grid->DeleteDataPtr("gradTemperatureFaceX");
        grid->DeleteDataPtr("gradTemperatureFaceY");
        grid->DeleteDataPtr("gradTemperatureFaceZ");

        RDouble4D *gradUVWTCellCenterX = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
        RDouble4D *gradUVWTCellCenterY = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
        RDouble4D *gradUVWTCellCenterZ = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));
        delete gradUVWTCellCenterX; gradUVWTCellCenterX = 0;
        delete gradUVWTCellCenterY; gradUVWTCellCenterY = 0;
        delete gradUVWTCellCenterZ; gradUVWTCellCenterZ = 0;
        grid->DeleteDataPtr("gradUVWTCellCenterX");
        grid->DeleteDataPtr("gradUVWTCellCenterY");
        grid->DeleteDataPtr("gradUVWTCellCenterZ");

        RDouble4D *rotNSgradValueX = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotNSgradValueX"));
        RDouble4D *rotNSgradValueY = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotNSgradValueY"));
        RDouble4D *rotNSgradValueZ = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotNSgradValueZ"));
        delete rotNSgradValueX; rotNSgradValueX = 0;
        delete rotNSgradValueY; rotNSgradValueY = 0;
        delete rotNSgradValueZ; rotNSgradValueZ = 0;
        grid->DeleteDataPtr("rotNSgradValueX");
        grid->DeleteDataPtr("rotNSgradValueY");
        grid->DeleteDataPtr("rotNSgradValueZ");

        RDouble4D *ql1 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1"));
        RDouble4D *qr1 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1"));
        RDouble4D *ql2 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2"));
        RDouble4D *qr2 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2"));
        RDouble4D *ql3 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3"));
        RDouble4D *qr3 = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3"));
        delete ql1; ql1 = 0;
        delete qr1; qr1 = 0;
        delete ql2; ql2 = 0;
        delete qr2; qr2 = 0;
        delete ql3; ql3 = 0;
        delete qr3; qr3 = 0;
        grid->DeleteDataPtr("ql1");
        grid->DeleteDataPtr("qr1");
        grid->DeleteDataPtr("ql2");
        grid->DeleteDataPtr("qr2");
        grid->DeleteDataPtr("ql3");
        grid->DeleteDataPtr("qr3");

        Int1D *speciesOrder = reinterpret_cast <Int1D *> (grid->GetDataPtr("speciesOrder"));
        delete speciesOrder; speciesOrder = 0;
        grid->DeleteDataPtr("speciesOrder");

        int iCodeBranch = GlobalDataBase::GetIntParaFromDB("iCodeBranch");
        if (iCodeBranch > 0)
        {
            RDouble *maxResidualVariation = reinterpret_cast <RDouble *> (grid->GetDataPtr("maxResidualVariation"));
            delete[]maxResidualVariation;
            grid->DeleteDataPtr("maxResidualVariation");
        }

        int nChemical = parameters->GetChemicalFlag();
        int nSpeciesNumber = parameters->GetNumberOfSpecies();
        if (nChemical > 0 && nSpeciesNumber > 0)
        {
            RDouble4D *gradSpeciesCellCenterX = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dcdx"));
            RDouble4D *gradSpeciesCellCenterY = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dcdy"));
            RDouble4D *gradSpeciesCellCenterZ = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dcdz"));
            delete gradSpeciesCellCenterX; gradSpeciesCellCenterX = 0;
            delete gradSpeciesCellCenterY; gradSpeciesCellCenterY = 0;
            delete gradSpeciesCellCenterZ; gradSpeciesCellCenterZ = 0;
            grid->DeleteDataPtr("dcdx");
            grid->DeleteDataPtr("dcdy");
            grid->DeleteDataPtr("dcdz");

            RDouble4D *speciesEnthalpy = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("hSpecies"));
            delete speciesEnthalpy; speciesEnthalpy = 0;
            grid->DeleteDataPtr("hSpecies");
            if (nWallBC > 0)
            {
                RDouble3D *surfaceMassFraction = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
                delete surfaceMassFraction; surfaceMassFraction = 0;
            }
            grid->DeleteDataPtr("surfaceMassFraction");

            RDouble3D *tflag = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("tflag"));
            delete tflag; tflag = 0;
            grid->DeleteDataPtr("tflag");

            RDouble3D *totalCvtr = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCvtr"));
            RDouble3D *totalCvv = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCvv"));
            RDouble3D *totalCve = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCve"));
            RDouble4D *speciesCvs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesCvs"));
            RDouble4D *speciesCps = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesCps"));
            RDouble4D *speciesCvvs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesCvvs"));
            RDouble4D *speciesCves = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesCves"));
            RDouble4D *speciesEtrs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEtrs"));
            RDouble4D *speciesEvs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEvs"));
            RDouble4D *speciesEes = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEes"));
            RDouble4D *speciesEs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEs"));
            delete totalCvtr; totalCvtr = 0;
            delete totalCvv; totalCvv = 0;
            delete totalCve; totalCve = 0;
            delete speciesCvs; speciesCvs = 0;
            delete speciesCps; speciesCps = 0;
            delete speciesCvvs; speciesCvvs = 0;
            delete speciesCves; speciesCves = 0;
            delete speciesEtrs; speciesEtrs = 0;
            delete speciesEvs; speciesEvs = 0;
            delete speciesEes; speciesEes = 0;
            delete speciesEs; speciesEs = 0;
            grid->DeleteDataPtr("totalCvtr");
            grid->DeleteDataPtr("totalCvv");
            grid->DeleteDataPtr("totalCve");
            grid->DeleteDataPtr("speciesCvs");
            grid->DeleteDataPtr("speciesCps");
            grid->DeleteDataPtr("speciesCvvs");
            grid->DeleteDataPtr("speciesCves");
            grid->DeleteDataPtr("speciesEtrs");
            grid->DeleteDataPtr("speciesEvs");
            grid->DeleteDataPtr("speciesEes");
            grid->DeleteDataPtr("speciesEs");

            RDouble3D *totalEnergy = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalEnergy"));
            RDouble3D *totalEnthalpy = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalEnthalpy"));
            RDouble3D *totalCp = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCp"));
            RDouble3D *totalCv = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCv"));
            delete totalEnergy; totalEnergy = 0;
            delete totalEnthalpy; totalEnthalpy = 0;
            delete totalCp; totalCp = 0;
            delete totalCv; totalCv = 0;
            grid->DeleteDataPtr("totalEnergy");
            grid->DeleteDataPtr("totalEnthalpy");
            grid->DeleteDataPtr("totalCp");
            grid->DeleteDataPtr("totalCv");

            int nTemperatureModel = parameters->GetTemperatureModel();
            int nEnergyRecycle = parameters->GetnEnergyRecycle();
            if (nTemperatureModel > 1 && nEnergyRecycle == 2)
            {
                RDouble4D *temporaryEve = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("temporaryEve"));
                delete temporaryEve; temporaryEve = 0;   //! temporaryEve.
                grid->DeleteDataPtr("temporaryEve");
            }
        }

        if (nWallBC > 0)
        {
            RDouble2D *firstLayerHeight = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("firstLayerHeight"));
            delete firstLayerHeight; firstLayerHeight = 0;

            RDouble2D *temperatureWall = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));
            delete temperatureWall; temperatureWall = 0;
            grid->DeleteDataPtr("surfaceTemperature");
        }
        grid->DeleteDataPtr("firstLayerHeight");

        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        if (nSlipBCModel > 0)
        {
            if (nWallBC > 0)
            {
                RDouble3D *surfaceSlipVariables = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceSlipVariables"));
                delete surfaceSlipVariables; surfaceSlipVariables = 0;
            }
            grid->DeleteDataPtr("surfaceSlipVariables");
        }
        int isAblationWall = GlobalDataBase::GetIntParaFromDB("nAblation");
        int Injection = GlobalDataBase::GetIntParaFromDB("isInjection");
        if (isAblationWall > 0 && Injection > 0)
        {
            if (nWallBC > 0)
            {
                RDouble2D *injectionVelocity = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("injectionVelocity"));
                delete injectionVelocity; injectionVelocity = 0;
                grid->DeleteDataPtr("injectionVelocity");
            }
        }

        int nRapidFlowfield = parameters->GetRapidFlowfieldMethod();
        if (nRapidFlowfield > 0)
        {
            if (nWallBC > 0)
            {
                RDouble3D *surfacePressure = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfacePressure"));
                delete surfacePressure; surfacePressure = 0;
            }
            grid->DeleteDataPtr("surfacePressure");
        }

        int nSurfHeatMonitor = parameters->GetSurfaceHeatingMonitor();
        if (nSurfHeatMonitor > 0)
        {
            if (nWallBC > 0)
            {
                RDouble3D *surfaceHeatFlux = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceHeatFlux"));
                delete surfaceHeatFlux; surfaceHeatFlux = 0;
            }
            grid->DeleteDataPtr("surfaceHeatFlux");
        }

        nWallBC = 1; //! The number of surface boundary conditions.
        nMaxSurfaceCell = 1; //! The number of surface cells.
        nIonizationFlag = 0; //! The flag of chemical reactions.
    }

    void NSSolverStruct::GetSurfaceCellNumber(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        //! Obtain the number of boundary condition regions.
        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        nWallBC = 0;
        nMaxSurfaceCell = 0;
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();

            if (!IsWall(BCType) && BCType != PHENGLEI::ABLATION_SURFACE)
            {
                continue;
            }
            ++ nWallBC; //The number of surface regions.

            //The following are surface boundary codes.
            int ist, ied, jst, jed, kst, ked;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            int ni = ied - ist + 1;
            int nj = jed - jst + 1;
            int nk = ked - kst + 1;
            int nCell = ni * nj * nk;
            nMaxSurfaceCell = MAX(nCell, nMaxSurfaceCell);
        }
    }

    void NSSolverStruct::InitializeSurfaceParameters(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        Param_NSSolverStruct *parameters = GetControlParameters();

        //! Set surface temperatures.
        RDouble Twall = parameters->GetWallTemperature();
        RDouble refTemperature = parameters->GetRefDimensionalTemperature();
        RDouble refT = Twall / refTemperature;
        if (nWallBC > 0)
        {
            if (Twall <= EPSILON)
            {
                refT = 1.0; //! Radiation equilibrium temperature wall.
            }

            RDouble2D *temperatureWall = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));
            if (temperatureWall == NULL)
            {
                temperatureWall = new RDouble2D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), fortranArray);
                grid->UpdateDataPtr("surfaceTemperature", temperatureWall);
            }

            //! Initialization with the freestream temperature.
            for (int iWall = 0; iWall < nWallBC; ++ iWall)
            {
                for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
                {
                    (*temperatureWall)(iWall, iCell) = refT;
                }
            }
        }

        //! Set surface mass fractions of species.
        RDouble initMassFraction[MAX_SPECIES_NUM] = { 0 };
        int nChemical = parameters->GetChemicalFlag();
        int nSpeciesNumber = parameters->GetNumberOfSpecies();
        int nNSEquation = parameters->GetNSEquationNumber();

        if (nChemical > 0 && nSpeciesNumber > 0)
        {
            int nIndexOfElectron = parameters->GetIndexOfElectron();
            if (nIndexOfElectron >= 0)
            {
                this->nIonizationFlag = 1;
            }

            RDouble *fcwMassFraction = parameters->GetFullyCatalyticMassFraction();
            for (int s = 0; s < nSpeciesNumber; ++ s)
            {
                initMassFraction[s] = fcwMassFraction[s];
            }

            RDouble catalyticCoef = parameters->GetCatalyticCoef();
            if (ABS(catalyticCoef) <= EPSILON) //! Initialization of the initial mass fraction of species.
            {
                RDouble *primitiveVar = parameters->GetPrimitiveVarFarfield();
                for (int s = 0; s < nSpeciesNumber; ++ s)
                {
                    initMassFraction[s] = primitiveVar[nNSEquation + s];
                }
            }

            RDouble3D *surfaceMassFraction = NULL;
            if (nWallBC > 0)
            {
                surfaceMassFraction = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
                if (surfaceMassFraction == NULL)
                {
                    surfaceMassFraction = new RDouble3D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), Range(0, nSpeciesNumber - 1), fortranArray);
                    grid->UpdateDataPtr("surfaceMassFraction", surfaceMassFraction);
                }
            }

            for (int iWall = 0; iWall < nWallBC; ++ iWall)
            {
                for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
                {
                    for (int s = 0; s < nSpeciesNumber; ++ s)
                    {
                        (*surfaceMassFraction)(iWall, iCell, s) = initMassFraction[s];
                    }
                }
            }
        }

        //! Set surface slip parameters.
        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        if (nSlipBCModel > 0)
        {
            //! The array includes variables such as slip temperatures(Tts, Tvs, Tes), slip velocity(Vx, Vy, Vz) and slip mass fractions of species.
            int nSlipVar = nSpeciesNumber + 6;
            RDouble3D *surfaceSlipVariables = NULL;
            if (nWallBC > 0)
            {
                surfaceSlipVariables = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceSlipVariables"));
                if (surfaceSlipVariables == NULL)
                {
                    surfaceSlipVariables = new RDouble3D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), Range(0, nSlipVar - 1), fortranArray);
                    grid->UpdateDataPtr("surfaceSlipVariables", surfaceSlipVariables);
                }
            }
            for (int iWall = 0; iWall < nWallBC; ++ iWall)
            {
                for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
                {
                    (*surfaceSlipVariables)(iWall, iCell, 0) = refT; //! Temperature.
                    (*surfaceSlipVariables)(iWall, iCell, 1) = refT;
                    (*surfaceSlipVariables)(iWall, iCell, 2) = refT;
                    (*surfaceSlipVariables)(iWall, iCell, 3) = 0.0; //! Velocity.
                    (*surfaceSlipVariables)(iWall, iCell, 4) = 0.0;
                    (*surfaceSlipVariables)(iWall, iCell, 5) = 0.0;
                    for (int s = 6; s < nSlipVar; ++ s) //! Species mass fractions.
                    {
                        (*surfaceSlipVariables)(iWall, iCell, s) = initMassFraction[s - 6];
                    }
                }
            }
        }
    }

    void NSSolverStruct::InitializeSurfaceParametersFromFile(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        Param_NSSolverStruct *parameters = GetControlParameters();
        RDouble refTemperature = parameters->GetRefDimensionalTemperature();

        //! Obtain the number of boundary condition regions.
        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        int nIndexOfWall = 0, nIndexOfCell = 0;

        RDouble2D *temperatureWall = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();

            if (!IsWall(BCType) && BCType != PHENGLEI::ABLATION_SURFACE)
            {
                continue;
            }

            RDouble3D *wallTempArray;
            if (structBC->CheckFieldData("wallTempArray"))
            {
                wallTempArray = reinterpret_cast <RDouble3D *> (structBC->GetFieldDataPtr("wallTempArray"));
            }
            else
            {
                continue;
            }

            //! The following are surface boundary codes.
            int ist, ied, jst, jed, kst, ked;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            nIndexOfCell = 0;
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*temperatureWall)(nIndexOfWall, nIndexOfCell) = (*wallTempArray)(i, j, k) / refTemperature;

                        ++ nIndexOfCell; //! Next grid cell.
                    }
                }
            }

            ++ nIndexOfWall; //! Next wall boundary.
        }
    }

    void NSSolverStruct::ChangeSurfaceTemperature(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        Param_NSSolverStruct *parameters = GetControlParameters();

        //! Set surface temperatures.
        RDouble Twall = parameters->GetWallTemperature();
        RDouble refTemperature = parameters->GetRefDimensionalTemperature();
        RDouble refT = Twall / refTemperature;
        if (nWallBC > 0)
        {
            if (Twall <= EPSILON)
            {
                refT = 1.0; //! Radiation equilibrium temperature wall.
            }

            RDouble2D *temperatureWall = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));

            //! Initialization with the freestream temperature.
            for (int iWall = 0; iWall < nWallBC; ++ iWall)
            {
                for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
                {
                    (*temperatureWall)(iWall, iCell) = refT;
                }
            }
        }
    }

    void NSSolverStruct::ChangeSurfaceMassFractions(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        Param_NSSolverStruct *parameters = GetControlParameters();

        //! Set surface mass fractions of species.
        RDouble initMassFraction[MAX_SPECIES_NUM];
        int nChemical = parameters->GetChemicalFlag();
        int nSpeciesNumber = parameters->GetNumberOfSpecies();
        int nNSEquation = parameters->GetNSEquationNumber();

        if (nChemical > 0 && nSpeciesNumber > 0)
        {
            RDouble *fcwMassFraction = parameters->GetFullyCatalyticMassFraction();
            for (int s = 0; s < nSpeciesNumber; ++ s)
            {
                initMassFraction[s] = fcwMassFraction[s];
            }

            RDouble catalyticCoef = parameters->GetCatalyticCoef();
            if (ABS(catalyticCoef) <= EPSILON) //! Initialization of the initial mass fraction of species.
            {
                RDouble *primitiveVar = parameters->GetPrimitiveVarFarfield();
                for (int s = 0; s < nSpeciesNumber; ++ s)
                {
                    initMassFraction[s] = primitiveVar[nNSEquation + s];
                }
            }

            RDouble3D *surfaceMassFraction = NULL;
            if (nWallBC > 0)
            {
                surfaceMassFraction = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
            }

            for (int iWall = 0; iWall < nWallBC; ++ iWall)
            {
                for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
                {
                    for (int s = 0; s < nSpeciesNumber; ++ s)
                    {
                        (*surfaceMassFraction)(iWall, iCell, s) = initMassFraction[s];
                    }
                }
            }
        }
    }

    bool NSSolverStruct::JudgeIfRestart()
    {
        string restartNSFile = ".\results\flow.dat";
        GlobalDataBase::GetData("restartNSFile", &restartNSFile, PHSTRING, 1);

        if (PHMPI::IsParallelRun())
        {
            restartNSFile = PHSPACE::AddSymbolToFileName(restartNSFile, "_", 0);
        }

        bool restart_flag = false;

        ifstream infile(restartNSFile.c_str(), ios::in);
        if (infile)
        {
            restart_flag = true;
            infile.close();
            infile.clear();
        }

        return restart_flag;
    }

    bool NSSolverStruct::JudgeIfReadAverage()
    {
        string restartNSFile = ".\results\flow.dat";
        GlobalDataBase::GetData("restartNSFile", &restartNSFile, PHSTRING, 1);
        string averageFlowFile = AddSymbolToFileName(restartNSFile, "_", "Average");

        if (PHMPI::IsParallelRun())
        {
            averageFlowFile = PHSPACE::AddSymbolToFileName(averageFlowFile, "_", 0);
        }

        bool restart_flag = false;

        ifstream infile(averageFlowFile.c_str(), ios::in);
        if (infile)
        {
            restart_flag = true;
            infile.close();
            infile.clear();
        }

        return restart_flag;
    }

    void NSSolverStruct::DumpRestart(ActionKey *actkey)
    {
        int level = actkey->level;
        StructGrid *grid = StructGridCast(GetGrid(level));

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nEquation = GetNumberOfEquations();

        Range M(0, nEquation - 1);
        int mst = M.first();
        int med = M.last();

        int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");

        DataContainer *restartData = actkey->GetData();
        restartData->moveto_begin();
        restartData->write(&outIterStep, sizeof(int));

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        restartData->write(&q(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        if (nChemical == 1)
        {
            for (int m = 0; m < nTemperatureModel; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            restartData->write(&temperatures(i, j, k, m), sizeof(RDouble));
                        }
                    }
                }
            }
        }

        int isUnsteady = parameters->GetIsUnsteady();
        if (!isUnsteady)
        {
            return;
        }

        RDouble physicalTime = GlobalDataBase::GetDoubleParaFromDB("physicalTime");
        restartData->write(&physicalTime, sizeof(RDouble));

        RDouble4D &qn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
        RDouble4D &qn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        restartData->write(&qn1(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        restartData->write(&qn2(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &residualn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n1"));
        RDouble4D &residualn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n2"));

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        restartData->write(&residual(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        restartData->write(&residualn1(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        restartData->write(&residualn2(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }
    }

    void NSSolverStruct::DumpRestartH5(ActionKey *actkey)
    {
        using namespace PHMPI;
        int currentProcessorID = GetCurrentProcessorID();

        int version = 1;
        Param_NSSolverStruct *parameters = GetControlParameters();
        int isUnsteady = parameters->GetIsUnsteady();
        int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
        if (currentProcessorID == GetServerProcessorID())
        {
            WriteData(actkey->filepos, &version, "Version");
            WriteData(actkey->filepos, &outIterStep, "outnstep");

            if (isUnsteady == 1)
            {
                RDouble physicalTime = GlobalDataBase::GetDoubleParaFromDB("physicalTime");
                WriteData(actkey->filepos, &physicalTime, "physicalTime");

                if (IsNeedStatistics())
                {
                    int nStatisticalStep = 0;
                    GlobalDataBase::GetData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
                    WriteData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
                }
            }
        }

        Grid *grid = GetGrid(actkey->level);
        int gridID = grid->GetZoneID();

        StructGrid *strGrid = StructGridCast(GetGrid(actkey->level));
        int nTotalCell = strGrid->GetNTotalCell();

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        int nEquation = GetNumberOfEquations();

        Range M(0, nEquation - 1);
        int mst = M.first();
        int med = M.last();

        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        //RFloat * primitiveVariables = &q(ist, jst, kst, mst);

        hid_t grploc;
        string grpName;

        ostringstream oss;
        oss << "Group" << gridID;
        grpName = oss.str();

        grploc = OpenGroup(actkey->filepos, grpName);

        //WriteData(grploc, primitiveVariables, "q");
        WriteData(grploc, &primitiveVariables(iCellStart, jCellStart, kCellStart, 0), "q");
        WriteData(grploc, &nTotalCell, "nTotalCell");

        //=======dump data for surface info=======
        RDouble wallTemperature = parameters->GetWallTemperature();
        RDouble catalyticCoef = parameters->GetCatalyticCoef();
        int nChemical = parameters->GetChemicalFlag();
        int numberOfSpecies = parameters->GetNumberOfSpecies();
        int nTemperatureModel = parameters->GetTemperatureModel();

        WriteData(grploc, &wallTemperature, "surfaceTemperatureFlag");
        WriteData(grploc, &catalyticCoef, "catalyticCoefFlag");
        WriteData(grploc, &nChemical, "nIsChemicalFlow");
        WriteData(grploc, &numberOfSpecies, "numberOfSpecies");
        WriteData(grploc, &nTemperatureModel, "nTemperatureModel");
        Int1D &speciesOrder = *reinterpret_cast <Int1D *> (grid->GetDataPtr("speciesOrder"));
        WriteData(grploc, &speciesOrder(0), "speciesNameList");

        //dump temperatures to file.
        //if (nChemical > 0)
        //{
        WriteData(grploc, &temperatures(iCellStart, jCellStart, kCellStart, mst), "t");
        //}

        if (nWallBC > 0)    //! The wall boundary condition exists.
        {
            //! Write the temperature on wall.
            RDouble2D &temperatureWall = *reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));
            WriteData(grploc, &temperatureWall(0, 0), "surfaceTemperature");

            //! Write the mass fractions of species on wall.
            if (nChemical > 0 && numberOfSpecies > 0)
            {
                RDouble3D &surfaceMassFraction = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
                WriteData(grploc, &surfaceMassFraction(0, 0, 0), "surfaceMassFraction");
            }
        }

        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        if (nSlipBCModel > 0 && nWallBC > 0)
        {
            RDouble3D &surfaceParameterSlip = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceSlipVariables"));
            WriteData(grploc, &surfaceParameterSlip(0, 0, 0), "surfaceSlipVariables");
        }

        if (isUnsteady)
        {
            RDouble4D &qn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
            WriteData(grploc, &qn1(iCellStart, jCellStart, kCellStart, mst), "q_unsteady_n1");

            RDouble4D &qn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));
            WriteData(grploc, &qn2(iCellStart, jCellStart, kCellStart, mst), "q_unsteady_n2");

            RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
            WriteData(grploc, &residual(iCellStart, jCellStart, kCellStart, mst), "res");

            RDouble4D &residualn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n1"));
            WriteData(grploc, &residualn1(iCellStart, jCellStart, kCellStart, mst), "res_unsteady_n1");

            RDouble4D &residualn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n2"));
            WriteData(grploc, &residualn2(iCellStart, jCellStart, kCellStart, mst), "res_unsteady_n2");

            RDouble4D &qAverage = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage"));
            WriteData(grploc, &qAverage(iCellStart, jCellStart, kCellStart, mst), "qAverage");

            RDouble4D &tauAverage = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("tauAverage"));
            WriteData(grploc, &tauAverage(iCellStart, jCellStart, kCellStart, 0), "tauAverage");

            RDouble4D &q2Average = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q2Average"));
            WriteData(grploc, &q2Average(iCellStart, jCellStart, kCellStart, 0), "q2Average");
        }

        H5Gclose(grploc);
    }

    void NSSolverStruct::ReadRestartFlowH5(ActionKey *actkey)
    {
        Grid *grid = GetGrid(actkey->level);
        int gridID = grid->GetZoneID();

        StructGrid *strGrid = StructGridCast(grid);
        int nTotalCell = strGrid->GetNTotalCell();

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nEquation = GetNumberOfEquations();
        Range M(0, nEquation - 1);
        int mst = M.first();
        int med = M.last();

        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble4D &heatTransferCoeff = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("heatTransferCoeff"));
        Int1D &speciesOrder = *reinterpret_cast <Int1D *> (grid->GetDataPtr("speciesOrder"));
        heatTransferCoeff = 0.0;

        int outIterStep = 0;
        ReadData(actkey->filepos, &outIterStep, "outnstep");
        GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

        hid_t grploc;
        string grpName;

        ostringstream oss;
        oss << "Group" << gridID;
        grpName = oss.str();

        grploc = OpenGroup(actkey->filepos, grpName);

        int nTotalCellRestart = 0;
        ReadData(grploc, &nTotalCellRestart, "nTotalCell");
        if (nTotalCellRestart != nTotalCell)
        {
            ostringstream erroeInfo;
            erroeInfo << " Error: the cell number in flow.dat is not equal to the cell number in grid file !" << endl;
            TK_Exit::ExceptionExit(erroeInfo.str());
        }

        //! add the data for surface info.
        int nChemical = parameters->GetChemicalFlag();
        int nNSEqn = parameters->GetNSEquationNumber();
        RDouble wallTemperature = parameters->GetWallTemperature();
        RDouble catalyticCoef = parameters->GetCatalyticCoef();
        int numberOfSpecies = parameters->GetNumberOfSpecies();
        int nTemperatureModel = parameters->GetTemperatureModel();
        RDouble preWallTemperature = wallTemperature, preCatalyticCoef = catalyticCoef;
        int preChemical = nChemical, preTemperatureModel = nTemperatureModel, preSpeciesNumber = numberOfSpecies;

        if (CheckDataExist(grploc, "surfaceTemperatureFlag"))
        {
            ReadData(grploc, &preWallTemperature, "surfaceTemperatureFlag");
        }
        if (CheckDataExist(grploc, "catalyticCoefFlag"))
        {
            ReadData(grploc, &preCatalyticCoef, "catalyticCoefFlag");
        }
        if (CheckDataExist(grploc, "nIsChemicalFlow"))
        {
            ReadData(grploc, &preChemical, "nIsChemicalFlow");
        }
        if (CheckDataExist(grploc, "numberOfSpecies"))
        {
            ReadData(grploc, &preSpeciesNumber, "numberOfSpecies");
        }
        if (CheckDataExist(grploc, "nTemperatureModel"))
        {
            ReadData(grploc, &preTemperatureModel, "nTemperatureModel");
        }
        int nTotalNumber = 16;
        Int1D *speciesOrder2 = new Int1D(Range(0, nTotalNumber - 1), fortranArray);
        //Initialization
        for (int n = 0; n < nTotalNumber; ++ n)
        {
            (*speciesOrder2)(n) = speciesOrder(n);
        }
        if (CheckDataExist(grploc, "speciesNameList"))
        {
            ReadData(grploc, &(*speciesOrder2)(0), "speciesNameList");
        }

        //read temperatures of flowfield.
        //if (nChemical > 0)
        //{
        if (CheckDataExist(grploc, "t"))
        {
            if (preTemperatureModel == nTemperatureModel)
            {
                ReadData(grploc, &temperatures(iCellStart, jCellStart, kCellStart, mst), "t");
            }
            else
            {
                RDouble4D *tmpTemperature = new RDouble4D(Range(iCellStart, iCellEnd), Range(jCellStart, jCellEnd), Range(kCellStart, kCellEnd), Range(0, preTemperatureModel - 1), fortranArray);
                ReadData(grploc, &(*tmpTemperature)(iCellStart, jCellStart, kCellStart, mst), "t");
                SetFlowfieldTemperature(grid, preChemical, preTemperatureModel, tmpTemperature);
                delete tmpTemperature;
            }
        }
        //}

        //Read the flowfield data.
        if (preChemical == nChemical && preTemperatureModel == nTemperatureModel && IsEqualGasModel(&speciesOrder, speciesOrder2))
        {
            ReadData(grploc, &primitiveVariables(iCellStart, jCellStart, kCellStart, mst), "q");
        }
        else
        {
            int nEqn = nNSEqn + preSpeciesNumber + preTemperatureModel - 1;
            RDouble4D *tmpVariables = new RDouble4D(Range(iCellStart, iCellEnd), Range(jCellStart, jCellEnd), Range(kCellStart, kCellEnd), Range(0, nEqn - 1), fortranArray);
            ReadData(grploc, &(*tmpVariables)(iCellStart, jCellStart, kCellStart, mst), "q");
            SetPrimitiveVariables(grid, preChemical, preTemperatureModel, preSpeciesNumber, speciesOrder2, tmpVariables);
            delete tmpVariables;
        }

        if (nWallBC > 0)    //! The wall boundary condition exists.
        {
            //! Write the temperature on wall.
            RDouble2D &temperatureWall = *reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));
            if (CheckDataExist(grploc, "surfaceTemperature"))
            {
                ReadData(grploc, &temperatureWall(0, 0), "surfaceTemperature");
                if (wallTemperature != preWallTemperature)
                {
                    ChangeSurfaceTemperature(grid);
                }
            }

            //! Write the mass fractions of species on wall.
            if (nChemical > 0 && numberOfSpecies > 0)
            {
                RDouble3D &surfaceMassFraction = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
                if (CheckDataExist(grploc, "surfaceMassFraction"))
                {
                    if (IsEqualGasModel(&speciesOrder, speciesOrder2))
                    {
                        ReadData(grploc, &surfaceMassFraction(0, 0, 0), "surfaceMassFraction");
                    }
                    else
                    {
                        RDouble3D *tmpSurfaceMassFraction = new RDouble3D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), Range(0, preSpeciesNumber - 1), fortranArray);
                        ReadData(grploc, &(*tmpSurfaceMassFraction)(0, 0, 0), "surfaceMassFraction");
                        SetSurfaceMassFractions(grid, preChemical, preSpeciesNumber, speciesOrder2, tmpSurfaceMassFraction);
                        delete tmpSurfaceMassFraction;
                    }
                }

                if (catalyticCoef != preCatalyticCoef)
                {
                    ChangeSurfaceMassFractions(grid);
                }
            }
        }

        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        if (nSlipBCModel > 0 && nWallBC > 0)
        {
            RDouble3D &surfaceParameterSlip = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceSlipVariables"));
            if (CheckDataExist(grploc, "surfaceSlipVariables"))
            {
                if (IsEqualGasModel(&speciesOrder, speciesOrder2))
                {
                    ReadData(grploc, &surfaceParameterSlip(0, 0, 0), "surfaceSlipVariables");
                }
                else
                {
                    int nLen = 6 + preSpeciesNumber;
                    RDouble3D *tmpSurfaceSlipVar = new RDouble3D(Range(0, nWallBC - 1), Range(0, nMaxSurfaceCell - 1), Range(0, nLen - 1), fortranArray);
                    ReadData(grploc, &(*tmpSurfaceSlipVar)(0, 0, 0), "surfaceSlipVariables");
                    SetSurfaceSlipVariables(grid, preChemical, preSpeciesNumber, speciesOrder2, tmpSurfaceSlipVar);
                    delete tmpSurfaceSlipVar;
                }
            }
        }

        //Remove dynamic memories.
        delete speciesOrder2;

        int isUnsteady = parameters->GetIsUnsteady();
        if (!isUnsteady)
        {
            RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            residual(i, j, k, m) = 0.0;
                        }
                    }
                }
            }

            H5Gclose(grploc);
            return;
        }

        RDouble4D &qn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
        RDouble4D &qn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));
        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &residualn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n1"));
        RDouble4D &residualn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n2"));
        RDouble4D &qAverage = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage"));
        RDouble4D &tauAverage = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("tauAverage"));
        RDouble4D &q2Average = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q2Average"));

        int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();
        RDouble physicalTime = 0.0;

        if (ifStartFromSteadyResults)
        {
            if (grid->GetGridID()->GetIndex() == 0)
            {
                PrintToWindow("Restart from steady NS flow field, reset outer step to be zero!\n");
            }

            //! Start from steady flow field.
            //! Reset the outer step when start from steady flow.
            outIterStep = 0;
            GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

            qn1 = primitiveVariables;
            qn2 = primitiveVariables;

            residual = 0.0;
            residualn1 = 0.0;
            residualn2 = 0.0;
        }
        else
        {
            ReadData(actkey->filepos, &physicalTime, "physicalTime");

            ReadData(grploc, &qn1(iCellStart, jCellStart, kCellStart, mst), "q_unsteady_n1");
            ReadData(grploc, &qn2(iCellStart, jCellStart, kCellStart, mst), "q_unsteady_n2");
            ReadData(grploc, &residual(iCellStart, jCellStart, kCellStart, mst), "res");
            ReadData(grploc, &residualn1(iCellStart, jCellStart, kCellStart, mst), "res_unsteady_n1");
            ReadData(grploc, &residualn2(iCellStart, jCellStart, kCellStart, mst), "res_unsteady_n2");
        }

        int nStatisticalStep = 0;

        if (IsNeedStatistics())
        {
            if (ifStartFromSteadyResults)
            {
                qAverage = primitiveVariables;
            }
            else
            {
                ReadData(grploc, &qAverage(iCellStart, jCellStart, kCellStart, mst), "qAverage");

                ReadData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
            }
        }

        if (IsNeedReynoldsStressStatistics())
        {
            if (ifStartFromSteadyResults)
            {
                tauAverage = 0.0;

                q2Average = 0.0;
            }
            else
            {
                ReadData(grploc, &tauAverage(iCellStart, jCellStart, kCellStart, 0), "tauAverage");

                ReadData(grploc, &q2Average(iCellStart, jCellStart, kCellStart, 0), "q2Average");
            }
        }

        GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

        GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);

        H5Gclose(grploc);
    }

    void NSSolverStruct::SetFlowfieldTemperature(Grid *grid, int nChemical, int nTemperatureModel, RDouble4D *primTemperature)
    {
        StructGrid *strGrid = StructGridCast(grid);
        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int curChem = parameters->GetChemicalFlag();
        int curTModel = parameters->GetTemperatureModel();

        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));

        if (curTModel <= nTemperatureModel) //The collection of the reading data is bigger.
        {
            for (int m = 0; m < curTModel; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            temperatures(i, j, k, m) = (*primTemperature)(i, j, k, m);
                        }
                    }
                }
            }
        }
        else //The collection of the reading data is smaller.
        {
            for (int m = 0; m < nTemperatureModel; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            temperatures(i, j, k, m) = (*primTemperature)(i, j, k, m);
                        }
                    }
                }
            }

            for (int m = nTemperatureModel; m < curTModel; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            temperatures(i, j, k, m) = (*primTemperature)(i, j, k, nTemperatureModel - 1);
                        }
                    }
                }
            }
        }
    }

    void NSSolverStruct::SetPrimitiveVariables(Grid *grid, int nChemical, int nTemperatureModel, int numberOfSpecies, Int1D *speciesOrder, RDouble4D *primVariables)
    {
        StructGrid *strGrid = StructGridCast(grid);
        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int curChem = parameters->GetChemicalFlag();
        int nNSEqn = parameters->GetNSEquationNumber();
        int curTModel = parameters->GetTemperatureModel();
        int curNumberOfSpecies = parameters->GetNumberOfSpecies();

        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
        RDouble *fcwMassFraction = parameters->GetFullyCatalyticMassFraction();
        Int1D &curSpeciesOrder = *reinterpret_cast <Int1D *> (grid->GetDataPtr("speciesOrder"));
        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble prim[MAX_SPECIES_NUM] = { 0.0 };

        //Set the variables of density, velocityU, velocityV, velocityW, and pressure.
        for (int m = 0; m < nNSEqn; ++ m)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        primitiveVariables(i, j, k, m) = (*primVariables)(i, j, k, m);
                    }
                }
            }
        }

        //Set the mass fractions of species.
        int nIndex1, nIndex2, totalNumber = 16;
        if (nChemical == 0) //the previous gas model is perfect gas.
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        for (int m = 0; m < curNumberOfSpecies; ++ m)
                        {
                            primitiveVariables(i, j, k, nNSEqn + m) = fcwMassFraction[m];
                        }
                    }
                }
            }
        }
        else //the previous gas model is real gas.
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        for (int m = 0; m < totalNumber; ++ m)
                        {
                            nIndex1 = curSpeciesOrder(m);
                            if (nIndex1 >= 0)
                            {
                                nIndex2 = (*speciesOrder)(m);
                                if (nIndex2 >= 0)
                                {
                                    prim[nNSEqn + nIndex1] = (*primVariables)(i, j, k, nNSEqn + nIndex2);
                                }
                                else
                                {
                                    //prim[nNSEqn + nIndex1] = primitiveVarFarfield[nNSEqn + nIndex1];
                                    prim[nNSEqn + nIndex1] = 0.0;
                                }
                            }
                        }

                        //Normalize.
                        if (curChem > 0)
                        {
                            gas->NormalizePrimitive(prim);
                            for (int m = 0; m < curNumberOfSpecies; ++ m)
                            {
                                primitiveVariables(i, j, k, nNSEqn + m) = prim[nNSEqn + m];
                            }
                        }
                    }
                }
            }
        }

        //Set the variables of vibrational energy and electron energy.
        if (curChem > 0 && curTModel > 1)
        {
            nIndex1 = nNSEqn + curNumberOfSpecies;
            nIndex2 = nNSEqn + numberOfSpecies;
            if (nChemical == 0) //the previous is perfect gas.
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            for (int m = 0; m < curTModel - 1; ++ m)
                            {
                                primitiveVariables(i, j, k, nIndex1 + m) = primitiveVarFarfield[nIndex1 + m];
                            }
                        }
                    }
                }
            }
            else //the previous is chemical flow.
            {
                if (curTModel <= nTemperatureModel) //The collection of the reading data is bigger.
                {
                    for (int m = 0; m < curTModel - 1; ++ m)
                    {
                        for (int k = kCellStart; k <= kCellEnd; ++ k)
                        {
                            for (int j = jCellStart; j <= jCellEnd; ++ j)
                            {
                                for (int i = iCellStart; i <= iCellEnd; ++ i)
                                {
                                    primitiveVariables(i, j, k, nIndex1 + m) = (*primVariables)(i, j, k, nIndex2 + m);
                                }
                            }
                        }
                    }
                }
                else //The collection of the reading data is smaller.
                {
                    for (int m = 0; m < nTemperatureModel - 1; ++ m)
                    {
                        for (int k = kCellStart; k <= kCellEnd; ++ k)
                        {
                            for (int j = jCellStart; j <= jCellEnd; ++ j)
                            {
                                for (int i = iCellStart; i <= iCellEnd; ++ i)
                                {
                                    primitiveVariables(i, j, k, nIndex1 + m) = (*primVariables)(i, j, k, nIndex2 + m);
                                }
                            }
                        }
                    }

                    RDouble vibEnergy = 0.0, eleEnergy = 0.0;
                    for (int k = kCellStart; k <= kCellEnd; ++ k)
                    {
                        for (int j = jCellStart; j <= jCellEnd; ++ j)
                        {
                            for (int i = iCellStart; i <= iCellEnd; ++ i)
                            {
                                for (int m = 0; m < curNumberOfSpecies; ++ m)
                                {
                                    prim[m] = primitiveVariables(i, j, k, nNSEqn + m);
                                }
                                RDouble temperature = temperatures(i, j, k, nTemperatureModel - 1);
                                if (temperature <= 0.0)
                                {
                                    temperature = 1.0;
                                }
                                vibEnergy = gas->GetMixedGasVibrationEnergy(temperature, prim);
                                eleEnergy = gas->GetMixedGasElectronEnergy(temperature, prim);

                                if (nTemperatureModel == 1)
                                {
                                    if (curTModel == 2)
                                    {
                                        primitiveVariables(i, j, k, nIndex1) = vibEnergy + eleEnergy;
                                    }
                                    else if (curTModel == 3)
                                    {
                                        primitiveVariables(i, j, k, nIndex1) = vibEnergy;
                                        primitiveVariables(i, j, k, nIndex1 + 1) = eleEnergy;
                                    }
                                }
                                else if (nTemperatureModel == 2)
                                {
                                    if (curTModel == 3)
                                    {
                                        primitiveVariables(i, j, k, nIndex1 + 1) = eleEnergy;
                                    }
                                }
                            }
                        }
                    }
                    //The collection of the reading data is bigger.
                }
                //the previous is chemical flow.
            }
            //Set vibrational energy and electron energy.
        }
    }

    void NSSolverStruct::SetSurfaceMassFractions(Grid *grid, int nChemical, int numberOfSpecies, Int1D *speciesOrder, RDouble3D *primMassFractions)
    {
        StructGrid *strGrid = StructGridCast(grid);
        Param_NSSolverStruct *parameters = GetControlParameters();
        int curChem = parameters->GetChemicalFlag();
        int nNSEqn = parameters->GetNSEquationNumber();
        int curNumberOfSpecies = parameters->GetNumberOfSpecies();
        RDouble catalytic = parameters->GetCatalyticCoef();

        RDouble *fcwMassFraction = parameters->GetFullyCatalyticMassFraction();
        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
        Int1D &curSpeciesOrder = *reinterpret_cast <Int1D *> (grid->GetDataPtr("speciesOrder"));
        RDouble3D &surfaceMassFraction = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
        RDouble prim[MAX_SPECIES_NUM] = { 0.0 };

        //Set the mass fractions of species.
        if (fabs(catalytic - 1.0) <= EPSILON) //fully catalytic wall condition.
        {
            for (int iWall = 0; iWall < nWallBC; ++ iWall)
            {
                for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
                {
                    for (int m = 0; m < curNumberOfSpecies; ++ m)
                    {
                        surfaceMassFraction(iWall, iCell, m) = fcwMassFraction[m];
                    }
                }
            }
        }
        else if (fabs(catalytic) <= EPSILON) //non-catalytic wall condition.
        {
            for (int iWall = 0; iWall < nWallBC; ++ iWall)
            {
                for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
                {
                    for (int m = 0; m < curNumberOfSpecies; ++ m)
                    {
                        surfaceMassFraction(iWall, iCell, m) = primitiveVarFarfield[nNSEqn + m];
                    }
                }
            }
        }
        else //partial catalytic wall condition.
        {
            int nIndex1, nIndex2, totalNumber = 16;
            for (int iWall = 0; iWall < nWallBC; ++ iWall)
            {
                for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
                {
                    for (int m = 0; m < totalNumber; ++ m)
                    {
                        nIndex1 = curSpeciesOrder(m);
                        if (nIndex1 >= 0)
                        {
                            nIndex2 = (*speciesOrder)(m);
                            if (nIndex2 >= 0)
                            {
                                prim[nNSEqn + nIndex1] = (*primMassFractions)(iWall, iCell, nIndex2);
                            }
                            else
                            {
                                prim[nNSEqn + nIndex1] = 0.0;
                            }
                        }
                    }

                    //Normalize.
                    gas->NormalizePrimitive(prim);
                    for (int m = 0; m < curNumberOfSpecies; ++ m)
                    {
                        surfaceMassFraction(iWall, iCell, m) = prim[nNSEqn + m];
                    }

                }
            }
            //partial catalytic wall condition.
        }
    }

    void NSSolverStruct::SetSurfaceSlipVariables(Grid *grid, int nChemical, int numberOfSpecies, Int1D *speciesOrder, RDouble3D *primSlipVar)
    {
        StructGrid *strGrid = StructGridCast(grid);
        Param_NSSolverStruct *parameters = GetControlParameters();
        int curChem = parameters->GetChemicalFlag();
        int nNSEqn = parameters->GetNSEquationNumber();
        int curNumberOfSpecies = parameters->GetNumberOfSpecies();

        Int1D &curSpeciesOrder = *reinterpret_cast <Int1D *> (grid->GetDataPtr("speciesOrder"));
        RDouble3D &surfaceParameterSlip = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceSlipVariables"));
        RDouble3D &surfaceMassFraction = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
        RDouble prim[MAX_SPECIES_NUM] = { 0.0 };

        //Set the temperature and velocity.
        int nIndex1, nIndex2, totalNumber = 6;
        for (int iWall = 0; iWall < nWallBC; ++ iWall)
        {
            for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
            {
                for (int m = 0; m < totalNumber; ++ m)
                {
                    surfaceParameterSlip(iWall, iCell, m) = (*primSlipVar)(iWall, iCell, m);
                }
            }
        }

        //Set the mass fractions.
        if (curChem > 0)
        {
            totalNumber = 16;
            for (int iWall = 0; iWall < nWallBC; ++ iWall)
            {
                for (int iCell = 0; iCell < nMaxSurfaceCell; ++ iCell)
                {
                    for (int m = 0; m < totalNumber; ++ m)
                    {
                        nIndex1 = curSpeciesOrder(m);
                        if (nIndex1 >= 0)
                        {
                            nIndex2 = (*speciesOrder)(m);
                            if (nIndex2 >= 0)
                            {
                                prim[nNSEqn + nIndex1] = (*primSlipVar)(iWall, iCell, nIndex2 + 6);
                            }
                            else
                            {
                                prim[nNSEqn + nIndex1] = surfaceMassFraction(iWall, iCell, nIndex1);
                            }
                        }
                    }

                    //Normalize.
                    gas->NormalizePrimitive(prim);
                    for (int m = 0; m < curNumberOfSpecies; ++ m)
                    {
                        surfaceParameterSlip(iWall, iCell, m + 6) = prim[nNSEqn + m];
                    }
                }
            }
        }
    }

    void NSSolverStruct::ReadPrefectGasRestartH5(ActionKey *actkey)
    {
        Grid *grid = GetGrid(actkey->level);
        int gridID = grid->GetZoneID();

        StructGrid *strGrid = StructGridCast(GetGrid(actkey->level));

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nEquation = GetNumberOfEquations();
        int nNSEquation = parameters->GetNSEquationNumber();
        Range M(0, nNSEquation - 1);
        int mst = M.first();
        int med = M.last();

        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        int outIterStep = 0;
        ReadData(actkey->filepos, &outIterStep, "outnstep");
        GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

        hid_t grploc;
        string grpName;

        ostringstream oss;
        oss << "Group" << gridID;
        grpName = oss.str();

        grploc = OpenGroup(actkey->filepos, grpName);

        ReadData(grploc, &primitiveVariables(iCellStart, jCellStart, kCellStart, mst), "q");

        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
        Range N(nNSEquation, nEquation - 1);
        int nst = N.first();
        int ned = N.last();

        for (int n = nst; n <= ned; ++ n)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        primitiveVariables(i, j, k, n) = primitiveVarFarfield[n];
                    }
                }
            }
        }

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        Range I, J, K;
        GetRange(I, J, K, iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd);
        residual(I, J, K, M) = 0.0;
        residual(I, J, K, N) = 0.0;

        H5Gclose(grploc);

        //! The unsteady flow is not considered here.
        //int isUnsteady = parameters->GetIsUnsteady();
    }

    void NSSolverStruct::ReadRestartH5(ActionKey *actkey)
    {
        int ifStartFromPerfectGasResults = GlobalDataBase::GetIntParaFromDB("ifStartFromPerfectGasResults");
        int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
        if (nChemical && ifStartFromPerfectGasResults)
        {
            ReadPrefectGasRestartH5(actkey);
            return;
        }

        int nContinueModel = 0;
        if (GlobalDataBase::IsExist("nContinueModel", PHINT, 1))
        {
            nContinueModel = GlobalDataBase::GetIntParaFromDB("nContinueModel");
        }
        if (nContinueModel == 1)
        {
            ReadRestartFlowH5(actkey);
            return;
        }

        Grid *grid = GetGrid(actkey->level);
        int gridID = grid->GetZoneID();

        StructGrid *strGrid = StructGridCast(GetGrid(actkey->level));
        int nTotalCell = strGrid->GetNTotalCell();

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nEquation = GetNumberOfEquations();
        Range M(0, nEquation - 1);
        int mst = M.first();
        int med = M.last();

        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble4D &heatTransferCoeff = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("heatTransferCoeff"));

        int outIterStep = 0;
        ReadData(actkey->filepos, &outIterStep, "outnstep");
        GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

        hid_t grploc;
        string grpName;

        ostringstream oss;
        oss << "Group" << gridID;
        grpName = oss.str();

        grploc = OpenGroup(actkey->filepos, grpName);

        int nTotalCellRestart = 0;
        ReadData(grploc, &nTotalCellRestart, "nTotalCell");
        if (nTotalCellRestart != nTotalCell)
        {
            ostringstream erroeInfo;
            erroeInfo << " Error: the cell number in flow.dat is not equal to the cell number in grid file !" << endl;
            TK_Exit::ExceptionExit(erroeInfo.str());
        }

        //Read the flowfield data.
        ReadData(grploc, &primitiveVariables(iCellStart, jCellStart, kCellStart, mst), "q");

        //! add the data for surface info.
        RDouble wallTemperature = parameters->GetWallTemperature();
        RDouble catalyticCoef = parameters->GetCatalyticCoef();
        int numberOfSpecies = parameters->GetNumberOfSpecies();
        int nTemperatureModel = parameters->GetTemperatureModel();
        RDouble preWallTemperature = wallTemperature, preCatalyticCoef = catalyticCoef;
        int preChemical = nChemical, preTemperatureModel = nTemperatureModel, preSpeciesNumber = numberOfSpecies;

        if (CheckDataExist(grploc, "surfaceTemperatureFlag"))
        {
            ReadData(grploc, &preWallTemperature, "surfaceTemperatureFlag");
        }
        if (CheckDataExist(grploc, "catalyticCoefFlag"))
        {
            ReadData(grploc, &preCatalyticCoef, "catalyticCoefFlag");
        }
        if (CheckDataExist(grploc, "nIsChemicalFlow"))
        {
            ReadData(grploc, &preChemical, "nIsChemicalFlow");
        }
        if (CheckDataExist(grploc, "numberOfSpecies"))
        {
            ReadData(grploc, &preSpeciesNumber, "numberOfSpecies");
        }
        if (CheckDataExist(grploc, "nTemperatureModel"))
        {
            ReadData(grploc, &preTemperatureModel, "nTemperatureModel");
        }

        //read temperatures.
        //if (nChemical > 0)
        //{
        if (CheckDataExist(grploc, "t"))
        {
            ReadData(grploc, &temperatures(iCellStart, jCellStart, kCellStart, mst), "t");
        }
        //}
        if (nWallBC > 0)    //! The wall boundary condition exists.
        {
            //! Write the temperature on wall.
            RDouble2D &temperatureWall = *reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));
            if (CheckDataExist(grploc, "surfaceTemperature"))
            {
                ReadData(grploc, &temperatureWall(0, 0), "surfaceTemperature");
                if (wallTemperature != preWallTemperature)
                {
                    ChangeSurfaceTemperature(grid);
                }
            }

            //! Write the mass fractions of species on wall.
            if (nChemical > 0 && numberOfSpecies > 0)
            {
                RDouble3D &surfaceMassFraction = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
                if (CheckDataExist(grploc, "surfaceMassFraction"))
                {
                    ReadData(grploc, &surfaceMassFraction(0, 0, 0), "surfaceMassFraction");
                }

                if (catalyticCoef != preCatalyticCoef)
                {
                    ChangeSurfaceMassFractions(grid);
                }
            }
        }

        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        if (nSlipBCModel > 0 && nWallBC > 0)
        {
            RDouble3D &surfaceParameterSlip = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceSlipVariables"));
            if (CheckDataExist(grploc, "surfaceSlipVariables"))
            {
                ReadData(grploc, &surfaceParameterSlip(0, 0, 0), "surfaceSlipVariables");
            }
        }

        int isUnsteady = parameters->GetIsUnsteady();
        if (!isUnsteady)
        {
            RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            residual(i, j, k, m) = 0.0;
                        }
                    }
                }
            }

            H5Gclose(grploc);
            return;
        }

        RDouble4D &qn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
        RDouble4D &qn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));
        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &residualn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n1"));
        RDouble4D &residualn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n2"));
        RDouble4D &qAverage = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage"));
        RDouble4D &tauAverage = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("tauAverage"));
        RDouble4D &q2Average = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q2Average"));

        int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();
        RDouble physicalTime = 0.0;

        if (ifStartFromSteadyResults)
        {
            if (grid->GetGridID()->GetIndex() == 0)
            {
                PrintToWindow("Restart from steady NS flow field, reset outer step to be zero!\n");
            }

            //! Start from steady flow field.
            //! Reset the outer step when start from steady flow.
            outIterStep = 0;
            GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

            qn1 = primitiveVariables;
            qn2 = primitiveVariables;

            residual = 0.0;
            residualn1 = 0.0;
            residualn2 = 0.0;
        }
        else
        {
            ReadData(actkey->filepos, &physicalTime, "physicalTime");

            ReadData(grploc, &qn1(iCellStart, jCellStart, kCellStart, mst), "q_unsteady_n1");
            ReadData(grploc, &qn2(iCellStart, jCellStart, kCellStart, mst), "q_unsteady_n2");
            ReadData(grploc, &residual(iCellStart, jCellStart, kCellStart, mst), "res");
            ReadData(grploc, &residualn1(iCellStart, jCellStart, kCellStart, mst), "res_unsteady_n1");
            ReadData(grploc, &residualn2(iCellStart, jCellStart, kCellStart, mst), "res_unsteady_n2");
        }

        int nStatisticalStep = 0;

        if (IsNeedStatistics())
        {
            if (ifStartFromSteadyResults)
            {
                qAverage = primitiveVariables;
            }
            else
            {
                ReadData(grploc, &qAverage(iCellStart, jCellStart, kCellStart, mst), "qAverage");

                ReadData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
            }
        }

        if (IsNeedReynoldsStressStatistics())
        {
            if (ifStartFromSteadyResults)
            {
                tauAverage = 0.0;

                q2Average = 0.0;
            }
            else
            {
                ReadData(grploc, &tauAverage(iCellStart, jCellStart, kCellStart, 0), "tauAverage");

                ReadData(grploc, &q2Average(iCellStart, jCellStart, kCellStart, 0), "q2Average");
            }
        }

        GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

        GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);

        H5Gclose(grploc);
    }

    void NSSolverStruct::ReadRestart(ActionKey *actkey)
    {
        int level = actkey->level;
        StructGrid *grid = StructGridCast(GetGrid(level));

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nEquation = GetNumberOfEquations();

        Range M(0, nEquation - 1);
        int mst = M.first();
        int med = M.last();

        DataContainer *restartData = actkey->GetData();
        restartData->moveto_begin();

        int outputStepInterval = 0;
        restartData->read(&outputStepInterval, sizeof(int));
        GlobalDataBase::UpdateData("outnstep", &outputStepInterval, PHINT, 1);

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        restartData->read(&q(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        if (nChemical == 1)
        {
            for (int m = 0; m < nTemperatureModel; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            restartData->read(&temperatures(i, j, k, m), sizeof(RDouble));
                        }
                    }
                }
            }
        }

        int isUnsteady = parameters->GetIsUnsteady();
        if (!isUnsteady)
        {
            RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            residual(i, j, k, m) = 0.0;
                        }
                    }
                }
            }
            return;
        }

        RDouble4D &qn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
        RDouble4D &qn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &residualn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n1"));
        RDouble4D &residualn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n2"));

        int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();

        RDouble physicalTime = 0.0;

        if (ifStartFromSteadyResults)
        {
            if (grid->GetGridID()->GetIndex() == 0)
            {
                PrintToWindow("Restart from steady NS flow field, reset outer step to be zero!\n");
            }

            //! Start from steady flow field.
            //! Reset the outer step when start from steady flow.
            outputStepInterval = 0;
            GlobalDataBase::UpdateData("outnstep", &outputStepInterval, PHINT, 1);

            qn1 = q;
            qn2 = q;

            residual = 0.0;
            residualn1 = 0.0;
            residualn2 = 0.0;
        }
        else
        {
            //! Start from unsteady flow field.
            restartData->read(&physicalTime, sizeof(RDouble));

            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            restartData->read(&qn1(i, j, k, m), sizeof(RDouble));
                        }
                    }
                }
            }

            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            restartData->read(&qn2(i, j, k, m), sizeof(RDouble));
                        }
                    }
                }
            }

            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            restartData->read(&residual(i, j, k, m), sizeof(RDouble));
                        }
                    }
                }
            }

            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            restartData->read(&residualn1(i, j, k, m), sizeof(RDouble));
                        }
                    }
                }
            }

            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            restartData->read(&residualn2(i, j, k, m), sizeof(RDouble));
                        }
                    }
                }
            }
        }

        GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);
    }

    void NSSolverStruct::AirForceCoef(ActionKey *actkey)
    {
        const int ALL_PART = -1;
        int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
        if (!isAle)
        {
            AirForceCoefParts(actkey, ALL_PART);

            int coordinate = LocalCoordinate;
            uint_t nTotalBC = GlobalBoundaryCondition::GetNumberOfBC();
            for (int iPart = 0; iPart < nTotalBC; ++ iPart)
            {
                if ((!GlobalBoundaryCondition::IsSolidWallBC(iPart)) && ((GlobalBoundaryCondition::GetBC(iPart)->GetBCType() != PHENGLEI::ABLATION_SURFACE)))
                {
                    continue;
                }

                //! Compute local part force.
                AirForceCoefParts(actkey, iPart, coordinate);

                //! Compute global part force.
                AirForceCoefParts(actkey, iPart);
            }
        }
        else
        {
            AirForceCoefComponents(actkey, ALL_PART);
        }
    }


    void NSSolverStruct::AirForceCoefParts(ActionKey *actkey, int partID, int coordinate)
    {
        int level = actkey->level;
        StructGrid *grid = StructGridCast(GetGrid(level));

        RDouble4D &primitiveVars = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));

        Param_NSSolverStruct *parameters = GetControlParameters();

        RDouble oRefReNumber = parameters->GetoRefReNumber();

        int viscousType = parameters->GetViscousType();

        const int ALL_PART = -1;
        SimpleBC *boundaryCondition = NULL;
		string partBCName = "";
        RDouble TorqueRefX = 0.0, TorqueRefY = 0.0, TorqueRefZ = 0.0;

        //! Compute hinge moment.
        int dumpHingeMoment = 0;
        RDouble localCoordAxis0[3] = { 0 };
        RDouble localCoordAxis1[3] = { 0 };

        if (partID == ALL_PART)
        {
            TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
            TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
            TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
        }
        else if ((partID != ALL_PART) && (coordinate == GlobalCoordinate))
        {
            boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
			partBCName = boundaryCondition->GetBCName();

            TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
            TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
            TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
        }
        else if ((partID != ALL_PART) && (coordinate == LocalCoordinate))
        {
            boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
			partBCName = boundaryCondition->GetBCName();
            boundaryCondition->GetParamData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
            boundaryCondition->GetParamData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
            boundaryCondition->GetParamData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);
       

            boundaryCondition->GetParamData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
            boundaryCondition->GetParamData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
            boundaryCondition->GetParamData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
        }

        using namespace IDX;
        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
        RDouble pressureFarfield = primitiveVarFarfield[IP];

        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());
        RDouble3D &vol = *(grid->GetCellVolume());

        RDouble cpx = 0.0;
        RDouble cpy = 0.0;
        RDouble cpz = 0.0;
        RDouble hingeMoment = 0.0;

        RDouble CA_f = 0.0;
        RDouble CA_p = 0.0;
        RDouble CN_f = 0.0;
        RDouble CN_p = 0.0;
        RDouble CZ_f = 0.0;
        RDouble CZ_p = 0.0;
        RDouble Cl_f = 0.0;
        RDouble Cl_p = 0.0;
        RDouble Cn_f = 0.0;
        RDouble Cn_p = 0.0;
        RDouble Cm_f = 0.0;
        RDouble Cm_p = 0.0;

        int ndim = GetDim();

        using namespace IDX;

        int ix = 0;

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *bcRegion = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = bcRegion->GetBCType();

            int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
            bcRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            bool isNeedComputeAirForce = (BCType == PHENGLEI::SOLID_SURFACE || BCType == PHENGLEI::ABLATION_SURFACE);
            if (isNeedComputeAirForce && partID != ALL_PART)
            {
                const string &bcName = bcRegion->GetBCName();
                isNeedComputeAirForce = isNeedComputeAirForce && (bcName == partBCName);
            }

            if (!isNeedComputeAirForce)
            {
                ix += (ied - ist + 1) * (jed - jst + 1) * (ked - kst + 1);
                continue;
            }

            int *s_lr3d = bcRegion->GetFaceDirectionIndex();
            int nSurface = bcRegion->GetFaceDirection() + 1;
            int leftOrRightIndex = bcRegion->GetFaceLeftOrRightIndex();

            int iWall = 0, jWall = 0, kWall = 0;

            RDouble dfx = 0, dfy = 0, dfz = 0, dpx = 0, dpy = 0, dpz = 0;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        //! i, j, k - is first cell index near SOLID_SURFACE.
                        bcRegion->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);

                        //! The surface direction of SOLID_SURFACE.
                        RDouble nx = leftOrRightIndex * xfv(iWall, jWall, kWall, nSurface);
                        RDouble ny = leftOrRightIndex * yfv(iWall, jWall, kWall, nSurface);
                        RDouble nz = leftOrRightIndex * zfv(iWall, jWall, kWall, nSurface);

                        //! The index of  first ghost cell layer near SOLID_SURFACE.
                        int iCellGhostLayer = i + s_lr3d[0];
                        int jCellGhostLayer = j + s_lr3d[1];
                        int kCellGhostLayer = k + s_lr3d[2];

                        RDouble pCellfirstLayer = primitiveVars(i, j, k, IP);

                        RDouble pressureWall = pCellfirstLayer;
                        RDouble pressureCoef = two * (pressureWall - pressureFarfield);

                        dpx = nx * pressureCoef;
                        dpy = ny * pressureCoef;
                        dpz = nz * pressureCoef;

                        CA_p += dpx;
                        CN_p += dpy;
                        CZ_p += dpz;

                        dfx = dpx;
                        dfy = dpy;
                        dfz = dpz;

                        RDouble xc = 0, yc = 0, zc = 0;
                        grid->FaceCoor(iWall, jWall, kWall, nSurface, xc, yc, zc);

                        RDouble visLaminar = 0, visTurbulence = 0;
                        RDouble dudx = 0, dudy = 0, dudz = 0, dvdx = 0, dvdy = 0, dvdz = 0, dwdx = 0, dwdy = 0, dwdz = 0;
                        if (viscousType > INVISCID)
                        {
                            Geometry3D *geometry3D = new Geometry3D;

                            if (nSurface == 1)
                            {
                                //! Get i-th face vector on the wall surface.
                                Get_GEO_I(i, j, k, iWall, iWall, ndim, geometry3D, xfv, yfv, zfv, vol);
                                //! Get derivative to the wall in i-direction of u,v,w.
                                DXDYDZ_I(primitiveVars, geometry3D, i, j, k, ndim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                DXDYDZ_I(primitiveVars, geometry3D, i, j, k, ndim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                DXDYDZ_I(primitiveVars, geometry3D, i, j, k, ndim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(iCellGhostLayer, j, k));
                                visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(iCellGhostLayer, j, k));
                            }
                            else if (nSurface == 2)
                            {
                                //! Get j-th face vector on the wall surface.
                                Get_GEO_J(i, j, k, jWall, jWall, ndim, geometry3D, xfv, yfv, zfv, vol);
                                //! Get derivative to the wall in j-direction of u,v,w.
                                DXDYDZ_J(primitiveVars, geometry3D, i, j, k, ndim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                DXDYDZ_J(primitiveVars, geometry3D, i, j, k, ndim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                DXDYDZ_J(primitiveVars, geometry3D, i, j, k, ndim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(i, jCellGhostLayer, k));
                                visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(i, jCellGhostLayer, k));
                            }
                            else
                            {
                                //! Get k-th face vector on the wall surface.
                                Get_GEO_K(i, j, k, kWall, kWall, ndim, geometry3D, xfv, yfv, zfv, vol);
                                //! Get derivative to the wall in k-direction of u,v,w.
                                DXDYDZ_K(primitiveVars, geometry3D, i, j, k, ndim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                DXDYDZ_K(primitiveVars, geometry3D, i, j, k, ndim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                DXDYDZ_K(primitiveVars, geometry3D, i, j, k, ndim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(i, j, kCellGhostLayer));
                                visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(i, j, kCellGhostLayer));
                            }

                            //! It is unnecessary to consider because the turbulence viscosity coefficient is 0.
                            //visTurbulence = 0.0;
                            RDouble viscousCoef = visLaminar + visTurbulence;

                            RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

                            //! Stress components.
                            RDouble txx = viscousCoef * (two * dudx - divv2p3);
                            RDouble tyy = viscousCoef * (two * dvdy - divv2p3);
                            RDouble tzz = viscousCoef * (two * dwdz - divv2p3);
                            RDouble txy = viscousCoef * (dudy + dvdx);
                            RDouble txz = viscousCoef * (dudz + dwdx);
                            RDouble tyz = viscousCoef * (dvdz + dwdy);

                            RDouble fvsx = -two * (nx * txx + ny * txy + nz * txz) * oRefReNumber;
                            RDouble fvsy = -two * (nx * txy + ny * tyy + nz * tyz) * oRefReNumber;
                            RDouble fvsz = -two * (nx * txz + ny * tyz + nz * tzz) * oRefReNumber;

                            dfx += fvsx;
                            dfy += fvsy;
                            dfz += fvsz;

                            CA_f += fvsx;
                            CN_f += fvsy;
                            CZ_f += fvsz;

                            Cl_f += (yc - TorqueRefY) * fvsz - (zc - TorqueRefZ) * fvsy;
                            Cn_f += (zc - TorqueRefZ) * fvsx - (xc - TorqueRefX) * fvsz;
                            Cm_f += (xc - TorqueRefX) * fvsy - (yc - TorqueRefY) * fvsx;

                            delete geometry3D;
                        }

                        Cl_p += (yc - TorqueRefY) * dpz - (zc - TorqueRefZ) * dpy;
                        Cn_p += (zc - TorqueRefZ) * dpx - (xc - TorqueRefX) * dpz;
                        Cm_p += (xc - TorqueRefX) * dpy - (yc - TorqueRefY) * dpx;

                        cpx += dpx;
                        cpy += dpy;
                        cpz += dpz;

                        if (dumpHingeMoment)
                        {
                            Post_ForceMoment forceMoment;
                            RDouble point[3] = { xc, yc, zc };
                            RDouble faceForce[3] = { dfx, dfy, dfz };
                            hingeMoment += forceMoment.ComputeMoment(point, localCoordAxis0, localCoordAxis1, faceForce);
                        }

                        ++ ix;
                    }
                }
            }
        }

        DataContainer *aerodynamicForceData = actkey->GetData();

        aerodynamicForceData->write(&CA_f, sizeof(RDouble));
        aerodynamicForceData->write(&CA_p, sizeof(RDouble));
        aerodynamicForceData->write(&CN_f, sizeof(RDouble));
        aerodynamicForceData->write(&CN_p, sizeof(RDouble));
        aerodynamicForceData->write(&CZ_f, sizeof(RDouble));
        aerodynamicForceData->write(&CZ_p, sizeof(RDouble));

        aerodynamicForceData->write(&cpx, sizeof(RDouble));
        aerodynamicForceData->write(&cpy, sizeof(RDouble));
        aerodynamicForceData->write(&cpz, sizeof(RDouble));

        aerodynamicForceData->write(&Cl_f, sizeof(RDouble));
        aerodynamicForceData->write(&Cl_p, sizeof(RDouble));
        aerodynamicForceData->write(&Cn_f, sizeof(RDouble));
        aerodynamicForceData->write(&Cn_p, sizeof(RDouble));
        aerodynamicForceData->write(&Cm_f, sizeof(RDouble));
        aerodynamicForceData->write(&Cm_p, sizeof(RDouble));

        aerodynamicForceData->write(&hingeMoment, sizeof(RDouble));
    }

    void NSSolverStruct::AirForceCoefComponents(ActionKey *actkey, int partID, int coordinate)
    {
        int level = actkey->level;
        StructGrid *grid = StructGridCast(GetGrid(level));

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &primitiveVars = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));

        RDouble4D &gradPrimtiveVarFaceX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
        RDouble4D &gradPrimtiveVarFaceY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
        RDouble4D &gradPrimtiveVarFaceZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

        RDouble4D &gradTemperatureFaceX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceX"));
        RDouble4D &gradTemperatureFaceY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceY"));
        RDouble4D &gradTemperatureFaceZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceZ"));

        RDouble4D &faceVelocityX = *reinterpret_cast <RDouble4D *> (grid->GetFaceVelocityX());
        RDouble4D &faceVelocityY = *reinterpret_cast <RDouble4D *> (grid->GetFaceVelocityY());
        RDouble4D &faceVelocityZ = *reinterpret_cast <RDouble4D *> (grid->GetFaceVelocityZ());

        Param_NSSolverStruct *parameters = GetControlParameters();

        RDouble oRefReNumber = parameters->GetoRefReNumber();

        int viscousType = parameters->GetViscousType();

        const int ALL_PART = -1;
        SimpleBC *boundaryCondition = 0;
        RDouble TorqueRefXGlobal = 0.0, TorqueRefYGlobal = 0.0, TorqueRefZGlobal = 0.0;
        RDouble TorqueRefXLocal = 0.0, TorqueRefYLocal = 0.0, TorqueRefZLocal = 0.0;

        TorqueRefXGlobal = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        TorqueRefYGlobal = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        TorqueRefZGlobal = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

        //! Compute hinge moment.
        int dumpHingeMoment = 0;
 //       RDouble localCoordAxis0[3];
 //       RDouble localCoordAxis1[3];

        using namespace IDX;
        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
        RDouble pressureFarfield = primitiveVarFarfield[IP];

        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());
        RDouble3D &vol = *(grid->GetCellVolume());
        RDouble4D &faceNormalVelocity = *(grid->GetFaceNormalVelocity());

        RDouble cpx = 0.0;
        RDouble cpy = 0.0;
        RDouble cpz = 0.0;
        RDouble hingeMoment = 0.0;

        int ndim = GetDim();

        using namespace IDX;

        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

        int numberOfAerodynamicForceComponents = GlobalDataBase::GetIntParaFromDB("numberOfAerodynamicForceComponents");
        vector< BasicAerodynamicForce * > *basicAerodynamicForceVector = CreateBasicAerodynamicForceVector();

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        int nTotalBC = GlobalBoundaryCondition::GetNumberOfBC();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *bcRegion = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = bcRegion->GetBCType();
            string bcName = bcRegion->GetBCName();

            int ist, ied, jst, jed, kst, ked;
            bcRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            bool isNeedComputeAirForce = (BCType == PHENGLEI::SOLID_SURFACE || BCType == PHENGLEI::ABLATION_SURFACE);
            if (!isNeedComputeAirForce)
            {
                continue;
            }

            int *faceDirection = bcRegion->GetFaceDirectionIndex();
            int nSurface = bcRegion->GetFaceDirection() + 1;
            int leftOrRightIndex = bcRegion->GetFaceLeftOrRightIndex();

            int iWall, jWall, kWall;
            int iCellGhostLayer, jCellGhostLayer, kCellGhostLayer;

            RDouble dfx, dfy, dfz, dpx, dpy, dpz;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        BasicAerodynamicForce aerodynamicForceGlobal;
                        for (int iPart = 0; iPart < nTotalBC; ++ iPart)
                        {
                            boundaryCondition = GlobalBoundaryCondition::GetBC(iPart);
                            string partName = boundaryCondition->GetBCName();
                            if (partName != bcName)
                            {
                                continue;
                            }
                            for (int iComponent = 0; iComponent < numberOfAerodynamicForceComponents; ++ iComponent)
                            {
                                if (FaceInAerodynamicForceComponents(bcName, iComponent))
                                {
                                    aerodynamicForceGlobal.SetMomentReferencePoint(TorqueRefXGlobal, TorqueRefYGlobal, TorqueRefZGlobal);
                                    //! i, j, k - is first face index near SOLID_SURFACE.
                                    bcRegion->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);

                                    //! The index of  first ghost cell layer near SOLID_SURFACE.
                                    bcRegion->GetGhostCellIndex(i, j, k, iCellGhostLayer, jCellGhostLayer, kCellGhostLayer, 1);

                                    //! The surface direction of SOLID_SURFACE.
                                    RDouble nx = leftOrRightIndex * xfv(iWall, jWall, kWall, nSurface);
                                    RDouble ny = leftOrRightIndex * yfv(iWall, jWall, kWall, nSurface);
                                    RDouble nz = leftOrRightIndex * zfv(iWall, jWall, kWall, nSurface);

                                    RDouble pCellfirstLayer = primitiveVars(i, j, k, IP);

                                    RDouble pressureWall = pCellfirstLayer;
                                    RDouble pressureCoef = two * (pressureWall - pressureFarfield);

                                    dpx = nx * pressureCoef;
                                    dpy = ny * pressureCoef;
                                    dpz = nz * pressureCoef;

                                    if (isUnsteady && isAle)
                                    {
                                        RDouble rhoWall = half * (primitiveVars(i, j, k, IR) + primitiveVars(iCellGhostLayer, jCellGhostLayer, kCellGhostLayer, IR));
                                        RDouble uWall = half * (primitiveVars(i, j, k, IU) + primitiveVars(iCellGhostLayer, jCellGhostLayer, kCellGhostLayer, IU));
                                        RDouble vWall = half * (primitiveVars(i, j, k, IV) + primitiveVars(iCellGhostLayer, jCellGhostLayer, kCellGhostLayer, IV));
                                        RDouble wWall = half * (primitiveVars(i, j, k, IW) + primitiveVars(iCellGhostLayer, jCellGhostLayer, kCellGhostLayer, IW));

                                        RDouble uGrid = faceVelocityX(iWall, jWall, kWall, nSurface);
                                        RDouble vGrid = faceVelocityY(iWall, jWall, kWall, nSurface);
                                        RDouble wGrid = faceVelocityZ(iWall, jWall, kWall, nSurface);

                                        RDouble gridVelocity = 2.0 * (uGrid * nx + vGrid * ny + wGrid * nz);

                                        dpx += rhoWall * uWall * gridVelocity;
                                        dpy += rhoWall * vWall * gridVelocity;
                                        dpz += rhoWall * wWall * gridVelocity;
                                    }

                                    aerodynamicForceGlobal.SetPressureAerodynamicForce(dpx, dpy, dpz);

                                    dfx = dpx;
                                    dfy = dpy;
                                    dfz = dpz;

                                    RDouble xc, yc, zc;
                                    grid->FaceCoor(iWall, jWall, kWall, nSurface, xc, yc, zc);

                                    aerodynamicForceGlobal.SetPressureAerodynamicMoment(xc, yc, zc);

                                    RDouble dudx = 0.0, dudy = 0.0, dudz = 0.0;
                                    RDouble dvdx = 0.0, dvdy = 0.0, dvdz = 0.0;
                                    RDouble dwdx = 0.0, dwdy = 0.0, dwdz = 0.0;
                                    RDouble visLaminar = 0.0, visTurbulence = 0.0;
                                    if (viscousType > INVISCID)
                                    {
                                        Geometry3D *geometry3D = new Geometry3D;

                                        if (nSurface == 1)
                                        {
                                            //! Get i-th face vector on the wall surface.
                                            Get_GEO_I(i, j, k, iWall, iWall, ndim, geometry3D, xfv, yfv, zfv, vol);
                                            //! Get derivative to the wall in i-direction of u,v,w.
                                            DXDYDZ_I(primitiveVars, geometry3D, i, j, k, ndim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                            DXDYDZ_I(primitiveVars, geometry3D, i, j, k, ndim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                            DXDYDZ_I(primitiveVars, geometry3D, i, j, k, ndim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                            visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(iCellGhostLayer, j, k));
                                            visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(iCellGhostLayer, j, k));
                                        }
                                        else if (nSurface == 2)
                                        {
                                            //! Get j-th face vector on the wall surface.
                                            Get_GEO_J(i, j, k, jWall, jWall, ndim, geometry3D, xfv, yfv, zfv, vol);
                                            //! Get derivative to the wall in j-direction of u,v,w.
                                            DXDYDZ_J(primitiveVars, geometry3D, i, j, k, ndim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                            DXDYDZ_J(primitiveVars, geometry3D, i, j, k, ndim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                            DXDYDZ_J(primitiveVars, geometry3D, i, j, k, ndim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                            visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(i, jCellGhostLayer, k));
                                            visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(i, jCellGhostLayer, k));
                                        }
                                        else
                                        {
                                            //! Get k-th face vector on the wall surface.
                                            Get_GEO_K(i, j, k, kWall, kWall, ndim, geometry3D, xfv, yfv, zfv, vol);
                                            //! Get derivative to the wall in k-direction of u,v,w.
                                            DXDYDZ_K(primitiveVars, geometry3D, i, j, k, ndim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                            DXDYDZ_K(primitiveVars, geometry3D, i, j, k, ndim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                            DXDYDZ_K(primitiveVars, geometry3D, i, j, k, ndim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                            visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(i, j, kCellGhostLayer));
                                            visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(i, j, kCellGhostLayer));
                                        }

                                        //! It is unnecessary to consider because the turbulence viscosity coefficient is 0.
                                        //visTurbulence = 0.0;
                                        RDouble viscousCoef = visLaminar + visTurbulence;

                                        RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

                                        //! Stress components.
                                        RDouble txx = viscousCoef * (two * dudx - divv2p3);
                                        RDouble tyy = viscousCoef * (two * dvdy - divv2p3);
                                        RDouble tzz = viscousCoef * (two * dwdz - divv2p3);
                                        RDouble txy = viscousCoef * (dudy + dvdx);
                                        RDouble txz = viscousCoef * (dudz + dwdx);
                                        RDouble tyz = viscousCoef * (dvdz + dwdy);

                                        RDouble fvsx = -two * (nx * txx + ny * txy + nz * txz) * oRefReNumber;
                                        RDouble fvsy = -two * (nx * txy + ny * tyy + nz * tyz) * oRefReNumber;
                                        RDouble fvsz = -two * (nx * txz + ny * tyz + nz * tzz) * oRefReNumber;

                                        aerodynamicForceGlobal.SetViscousAerodynamicForce(fvsx, fvsy, fvsz);

                                        aerodynamicForceGlobal.SetViscousAerodynamicMoment(xc, yc, zc);

                                        delete geometry3D;
                                    }
                                    aerodynamicForceGlobal.ComputeResultantAerodynamicForce();

                                    aerodynamicForceGlobal.ComputeResultantAerodynamicMoment();

                                    (*basicAerodynamicForceVector)[iComponent]->AddAerodynamicForce(&aerodynamicForceGlobal);

                                    (*basicAerodynamicForceVector)[numberOfAerodynamicForceComponents]->AddAerodynamicForce(&aerodynamicForceGlobal);
                                }
                            }
                        }
                    }
                }
            }
        }

        PHSPACE::WriteAerodynamicForceToActionKey(actkey, basicAerodynamicForceVector);
        PHSPACE::FreeBasicAerodynamicForceVector(basicAerodynamicForceVector);
    }

    //! Write the variables of Cp, Cf, Qw, y+ etc. on surface to the file.
    void NSSolverStruct::CpDistriCoef(ActionKey *actkey)
    {
        int level = actkey->level;
        DataContainer *wallDistributeData = actkey->GetData();
        wallDistributeData->moveto_begin();

        StructGrid *grid = StructGridCast(GetGrid(level));
        int GridID = grid->GetGridID()->GetIndex();
        PHWrite(wallDistributeData, GridID);

        if (level != 0) return;

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &primitiveVars = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperature = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));

        Param_NSSolverStruct *parameters = GetControlParameters();

        RDouble refReNumber = parameters->GetRefReNumber();
        RDouble oRefReNumber = parameters->GetoRefReNumber();

        RDouble refGama = parameters->GetRefGama();

        RDouble oPrandtlLaminar = parameters->GetoPrandtlLaminar();
        RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

        RDouble refMachNumber = parameters->GetRefMachNumber();
        int viscousType = parameters->GetViscousType();

        int mTT = parameters->GetmTT();
        int mTV = parameters->GetmTV();
        int mTE = parameters->GetmTE();

        using namespace IDX;
        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
        RDouble pressureFarfield = primitiveVarFarfield[IP];

        RDouble4D &gradPrimtiveVarFaceX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
        RDouble4D &gradPrimtiveVarFaceY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
        RDouble4D &gradPrimtiveVarFaceZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

        RDouble4D &gradTemperatureFaceX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceX"));
        RDouble4D &gradTemperatureFaceY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceY"));
        RDouble4D &gradTemperatureFaceZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceZ"));

        RDouble4D &heatTransferCoeff = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("heatTransferCoeff"));

        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());

        RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
        RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
        RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());
        RDouble4D &faceArea = *(grid->GetFaceArea());
        RDouble3D &cellVolume = *(grid->GetCellVolume());

        RDouble3D &xx = *grid->GetStructX();
        RDouble3D &yy = *grid->GetStructY();
        RDouble3D &zz = *grid->GetStructZ();

        int nsolid_surface = grid->GetNumberOfWallCell();

        if (nsolid_surface == 0)
        {
            return;
        }

        RDouble AoA = 0.0, attack = 0.0;
        AoA = parameters->GetAoA();
        attack = AoA * PI / 180.0;

        RDouble sina = sin(attack);
        RDouble cosa = cos(attack);

        int nDim = GetDim();
        int js = 1;
        int ks = 1;
        if (nDim == TWO_D)
        {
            ks = 0;
        }

        if (nDim == ONE_D)
        {
            js = 0;
        }

        RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);

        using namespace IDX;

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        //! Obtain the number of boundary condition regions.
        int nBCRegion = compositeBCRegion->GetnBCRegion();

        int coorCount, dataCount;

        int TecioMission = 0;

        //! Non-equilibrium Gas
        //! The universal gas constant.
        RDouble gasConstant = gas->GetUniversalGasConstant();
        //! The type of catalytic wall condition.
        RDouble catalyticCoef = parameters->GetCatalyticCoef();
        RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();
        int nChemical = parameters->GetChemicalFlag();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nTemperatureModel = parameters->GetTemperatureModel();
        RDouble Twall = parameters->GetWallTemperature();    //! Temperature on wall.
        RDouble referenceTemperature = parameters->GetRefDimensionalTemperature();
        RDouble2D *surfaceTemperature;
        if (Twall >= 0.0)
        {
            surfaceTemperature = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));
        }
        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        RDouble3D *surfaceSlipVariables = NULL;
        if (nSlipBCModel > 0)
        {
            surfaceSlipVariables = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceSlipVariables"));
        }
        RDouble2D *firstLayerHeight = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("firstLayerHeight"));
        int wallFunctionType = GlobalDataBase::GetIntParaFromDB("wallFunctionType");

        int numberOfSpecies = parameters->GetNumberOfSpecies();
        //! Allocate the dynamic memory space to save the primitive variables on wall.
        int nEquation = numberOfSpecies + nNSEquation + nTemperatureModel - 1;
        RDouble *wallVariables = new RDouble[nEquation];
        RDouble nondimTwall = Twall / referenceTemperature;
        RDouble wallTemperature[3] = { nondimTwall, nondimTwall, nondimTwall };
        RDouble Qw = 0.0, Qtr = 0.0, Qs = 0.0, Qv = 0.0, Qe = 0.0, Pw = 0.0, Tw = 0.0, Rhow = 0.0, absVelocity = 0.0;
        RDouble StantonNumber = 0.0, HeatCoef = 0.0, deltaT = 0.0, wallTotalEnthalpy = 0.0;

        int nSurfGradMethod = parameters->GetSurfaceGradientMethod();
        RDouble gradTx[3] = {0.0, 0.0, 0.0}, gradTy[3] = {0.0, 0.0, 0.0}, gradTz[3] = {0.0, 0.0, 0.0};
        RDouble gradCx[MAX_SPECIES_NUM] = {0.0}, gradCy[MAX_SPECIES_NUM] = {0.0}, gradCz[MAX_SPECIES_NUM] = {0.0};

        //! The mass fraction of oxygen and nitrogen under the fully catalytic condition.
        RDouble speciesEnthalpy[MAX_SPECIES_NUM], densityDiffusivity[MAX_SPECIES_NUM];
        string *speciesName;
        RDouble3D *surfaceMassFraction;
        if (nChemical > 0 && numberOfSpecies > 0)
        {
            speciesName = gas->GetNameOfSpecies();
            surfaceMassFraction = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
        }

        int nRapidFlowfield = parameters->GetRapidFlowfieldMethod();
        RDouble3D *surfacePressure = NULL;
        if (nRapidFlowfield > 0)
        {
            surfacePressure = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfacePressure"));
        }

        //! To obtain the reference values.
        RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
        RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
        RDouble refVelocity = refDimensionalSonicSpeed * refMachNumber;
        RDouble refSquareVelocity = refVelocity * refVelocity;
        RDouble refDynamicPressure = refDimensionalDensity * refSquareVelocity;
        RDouble refEnergy = refDimensionalDensity * refVelocity * refSquareVelocity;
        RDouble refHeatFluxDimension = refEnergy * 0.001; // The dimensional value of total heat flux (kW/m2).
        RDouble refTotalEnthalpy = gas->ComputeReferenceTotalEnthalpy();

        int *visualVariablesType = postVisualWall->GetVisualVariablesType();
        int nWallVariables = postVisualWall->GetVisualVariablesNumber();
        int dataNumber = nWallVariables;
        for (int m = 0; m < nWallVariables; ++ m)
        {
            if (visualVariablesType[m] == VISUAL_WALL_NS)
            {
                dataNumber += numberOfSpecies - 1;
            }
        }

        Range M(0, dataNumber - 1);
        RDouble visualVariables[50]; //Store the values of variables.
        //the order of variables in the computational array, eg. 0 denotes the coefficient of pressure,
        //and the value of array indicates the index of the variable in the computation array.
        int variablesOrder[50];
        for (int m = 0; m < 50; ++ m)
        {
            variablesOrder[m] = m;
        }

        Range I, J, K;
        GetRange(ni, nj, nk, -2, 1, I, J, K);

        int isWennScheme = parameters->GetWennSchemeFlag();
        if (isWennScheme == 1)
        {
            GetRange(ni, nj, nk, -3, 2, I, J, K);
        }

        RDouble4D dataStore(I, J, K, M, fortranArray);
        RDouble3D isWall(I, J, K, fortranArray);
        dataStore = 0.0;
        isWall = 0.0;

        int indexOfWall = 0, indexOfCell = 0;
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();

            if (!IsWall(BCType) && BCType != PHENGLEI::ABLATION_SURFACE)
            {
                continue;
            }

            int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            //! Obtain the surface information.
            int *s_lr3d = structBC->GetFaceDirectionIndex();
            int nSurface = structBC->GetFaceDirection() + 1;
            int leftOrRightIndex = structBC->GetFaceLeftOrRightIndex();

            indexOfCell = 0;
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        //! i, j, k - is first CELL index near SOLID_SURFACE.
                        //! iWall, jWall, kWall - is NODE index on SOLID_SURFACE.
                        int iWall, jWall, kWall;
                        structBC->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);

                        RDouble nx = leftOrRightIndex * xfv(iWall, jWall, kWall, nSurface);
                        RDouble ny = leftOrRightIndex * yfv(iWall, jWall, kWall, nSurface);
                        RDouble nz = leftOrRightIndex * zfv(iWall, jWall, kWall, nSurface);

                        ////! Obtain the centroid of the surface.
                        //RDouble xCenterWall, yCenterWall, zCenterWall;
                        //grid->FaceCoor(iWall, jWall, kWall, nSurface, xCenterWall, yCenterWall, zCenterWall);
                        ////! Obtain the centroid of the cell on the positive first layer.
                        //RDouble xCellCenter, yCellCenter, zCellCenter;
                        //grid->CenterCoor(i, j, k, xCellCenter, yCellCenter, zCellCenter);

                        RDouble normalComponetX = leftOrRightIndex * faceNormalComponentX(iWall, jWall, kWall, nSurface);
                        RDouble normalComponetY = leftOrRightIndex * faceNormalComponentY(iWall, jWall, kWall, nSurface);
                        RDouble normalComponetZ = leftOrRightIndex * faceNormalComponentZ(iWall, jWall, kWall, nSurface);
                        RDouble surfaceArea = faceArea(iWall, jWall, kWall, nSurface);

                        //! The projection of the dy called the distance between the center of the cell and that of the surface on normal vector.
                        //RDouble deltaHeight = fabs((xCellCenter - xCenterWall) * normalComponetX + (yCellCenter - yCenterWall) * normalComponetY + (zCellCenter - zCenterWall) * normalComponetZ);
                        RDouble deltaHeight = (*firstLayerHeight)(indexOfWall, indexOfCell);

                        //! The index of  first ghost cell layer near SOLID_SURFACE.
                        int iCellGhostLayer = i + s_lr3d[0];
                        int jCellGhostLayer = j + s_lr3d[1];
                        int kCellGhostLayer = k + s_lr3d[2];

                        RDouble pressureFirstCell = primitiveVars(i, j, k, IP);
                        wallVariables[IP] = pressureFirstCell;

                        //! To obtain the variables on wall.
                        if (nSlipBCModel > 0)
                        {
                            for (int m = 0; m < 3; ++ m)
                            {
                                wallTemperature[m] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m);
                                wallVariables[m + 1] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m + 3);
                            }
                            deltaT = wallTemperature[ITT] - nondimTwall;
                        }
                        else
                        {
                            if (viscousType > INVISCID)
                            {
                                wallVariables[IU] = 0.0;
                                wallVariables[IV] = 0.0;
                                wallVariables[IW] = 0.0;
                            }
                            else
                            {
                                //The velocity is used to extract the wall streamline, so this velocity is not the real velocity on wall.
                                wallVariables[IU] = primitiveVars(i, j, k, IU);
                                wallVariables[IV] = primitiveVars(i, j, k, IV);
                                wallVariables[IW] = primitiveVars(i, j, k, IW);
                            }

                            deltaT = 0.0;
                            //! Set temperature on wall.
                            if (Twall < 0.0) //Adabatic wall.
                            {
                                wallTemperature[ITT] = temperature(i, j, k, mTT);
                                wallTemperature[ITV] = temperature(i, j, k, mTV);
                                wallTemperature[ITE] = temperature(i, j, k, mTE);

                            }
                            else //Isothermal wall or radiation equilibrium temperature wall.
                            {
                                wallTemperature[ITT] = (*surfaceTemperature)(indexOfWall, indexOfCell);
                                wallTemperature[ITV] = wallTemperature[ITT];
                                wallTemperature[ITE] = wallTemperature[ITT];
                            }
                        }

                        //! Compute the density on wall using the state equation.
                        wallVariables[IR] = refGama * refMachNumber * refMachNumber * wallVariables[IP] / wallTemperature[ITT];
                        if (nChemical >= 1)
                        {
                            //! To obtain the mass fraction on wall.
                            if (ABS(catalyticCoef) <= EPSILON && BCType != PHENGLEI::ABLATION_SURFACE) // Fully non-catalytic wall.
                            {
                                for (int s = nNSEquation; s < numberOfSpecies + nNSEquation; ++ s)
                                {
                                    wallVariables[s] = primitiveVars(i, j, k, s);
                                }
                            }
                            else // Fully catalytic wall or finite catalytic wall.
                            {
                                for (int s = nNSEquation; s < numberOfSpecies + nNSEquation; ++ s)
                                {
                                    wallVariables[s] = (*surfaceMassFraction)(indexOfWall, indexOfCell, s - nNSEquation);
                                }
                            }

                            RDouble massReciprocal = 0.0;
                            if (nTemperatureModel > 1)     //! Multi-Temperature model.
                            {
                                RDouble ceDivideMe = 0.0;
                                massReciprocal = gas->ComputeMolecularWeightReciprocalWithoutElectron(wallVariables, ceDivideMe);
                                wallVariables[IR] = wallVariables[IP] / (gasConstant * wallTemperature[ITT] * massReciprocal + gasConstant * wallTemperature[ITE] * ceDivideMe);
                            }
                            else         //! One-Temperature model.
                            {
                                //! Obtain molecular weight of the mixture gas.
                                massReciprocal = gas->ComputeMolecularWeightReciprocal(wallVariables);
                                wallVariables[IR] = wallVariables[IP] / (gasConstant * wallTemperature[ITT] * massReciprocal);
                            }
                        }

                        //! The pressure on wall.
                        Rhow = wallVariables[IR] * refDimensionalDensity;
                        Pw = wallVariables[IP] * refDynamicPressure;
                        Tw = wallTemperature[ITT] * referenceTemperature;
                        absVelocity = wallVariables[IU] * wallVariables[IU] + wallVariables[IV] * wallVariables[IV] + wallVariables[IW] * wallVariables[IW];

                        wallTotalEnthalpy = 0.0;
                        gas->ComputeEnthalpyByPrimitive(wallVariables, refGama, wallTotalEnthalpy, wallTemperature);
                        wallTotalEnthalpy += 0.5 * absVelocity;
                        absVelocity = sqrt(absVelocity) * refVelocity;
                        wallTotalEnthalpy = wallTotalEnthalpy * refSquareVelocity;

                        //! Pressure drag.
                        RDouble pressureWall = pressureFirstCell;
                        RDouble pressCoef = two * (pressureWall - pressureFarfield);

                        RDouble dpx = nx * pressCoef;
                        RDouble dpy = ny * pressCoef;
                        RDouble dpz = nz * pressCoef;

                        RDouble dfx = dpx;
                        RDouble dfy = dpy;
                        RDouble dfz = dpz;

                        RDouble fvsx = 0.0;
                        RDouble fvsy = 0.0;
                        RDouble fvsz = 0.0;

                        RDouble heatFlux = zero;
                        RDouble viscousCoef;

                        RDouble dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
                        RDouble dtdx, dtdy, dtdz;
                        RDouble visLaminar = 0.0, visTurbulence = 0.0;
                        RDouble heatCoeff = 0.0;
                        if (viscousType > INVISCID)
                        {
                            Geometry3D *geometry3D = new Geometry3D;

                            if (nSurface == 1)
                            {
                                //! Get i-th face vector on the wall surface.
                                Get_GEO_I(i, j, k, iWall, iWall, nDim, geometry3D, xfv, yfv, zfv, cellVolume);
                                //! Get derivative to the wall in i-direction of u,v,w,T.
                                DXDYDZ_I(primitiveVars, geometry3D, i, j, k, nDim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                DXDYDZ_I(primitiveVars, geometry3D, i, j, k, nDim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                DXDYDZ_I(primitiveVars, geometry3D, i, j, k, nDim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                //DXDYDZ_I(temperature, geometry3D, i, j, k, nDim, ITT, dtdx, dtdy, dtdz, leftOrRightIndex);
                                for (int m = 0; m < nTemperatureModel; ++ m)
                                {
                                    DXDYDZ_I(temperature, geometry3D, i, j, k, nDim, m, gradTx[m], gradTy[m], gradTz[m], leftOrRightIndex);
                                }
                                for (int m = 0; m < numberOfSpecies; ++ m)
                                {
                                    DXDYDZ_I(primitiveVars, geometry3D, i, j, k, nDim, m + nNSEquation, gradCx[m], gradCy[m], gradCz[m], leftOrRightIndex);
                                }

                                visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(iCellGhostLayer, j, k));
                                visTurbulence = 0.0;
                                if (wallFunctionType != WALLFUNCTION::NONE)
                                {
                                    visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(iCellGhostLayer, j, k));
                                    heatCoeff = heatTransferCoeff(i, j, k, ITT) + heatTransferCoeff(iCellGhostLayer, j, k, ITT);
                                }
                            }
                            else if (nSurface == 2)
                            {
                                //! Get j-th face vector on the wall surface;.
                                Get_GEO_J(i, j, k, jWall, jWall, nDim, geometry3D, xfv, yfv, zfv, cellVolume);
                                //! Get derivative to the wall in j-direction of u,v,w,T.
                                DXDYDZ_J(primitiveVars, geometry3D, i, j, k, nDim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                DXDYDZ_J(primitiveVars, geometry3D, i, j, k, nDim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                DXDYDZ_J(primitiveVars, geometry3D, i, j, k, nDim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                //DXDYDZ_J(temperature, geometry3D, i, j, k, nDim, ITT, dtdx, dtdy, dtdz, leftOrRightIndex);
                                for (int m = 0; m < nTemperatureModel; ++ m)
                                {
                                    DXDYDZ_J(temperature, geometry3D, i, j, k, nDim, m, gradTx[m], gradTy[m], gradTz[m], leftOrRightIndex);
                                }
                                for (int m = 0; m < numberOfSpecies; ++ m)
                                {
                                    DXDYDZ_J(primitiveVars, geometry3D, i, j, k, nDim, m + nNSEquation, gradCx[m], gradCy[m], gradCz[m], leftOrRightIndex);
                                }

                                visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(i, jCellGhostLayer, k));
                                visTurbulence = 0.0;
                                if (wallFunctionType != WALLFUNCTION::NONE)
                                {
                                    visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(i, jCellGhostLayer, k));
                                    heatCoeff = heatTransferCoeff(i, j, k, ITT) + heatTransferCoeff(i, jCellGhostLayer, k, ITT);
                                }
                            }
                            else
                            {
                                //! Get k-th face vector on the wall surface.
                                Get_GEO_K(i, j, k, kWall, kWall, nDim, geometry3D, xfv, yfv, zfv, cellVolume);
                                //! Get derivative to the wall in k-direction of u,v,w,T.
                                DXDYDZ_K(primitiveVars, geometry3D, i, j, k, nDim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                DXDYDZ_K(primitiveVars, geometry3D, i, j, k, nDim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                DXDYDZ_K(primitiveVars, geometry3D, i, j, k, nDim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                //DXDYDZ_K(temperature, geometry3D, i, j, k, nDim, ITT, dtdx, dtdy, dtdz, leftOrRightIndex);
                                for (int m = 0; m < nTemperatureModel; ++ m)
                                {
                                    DXDYDZ_K(temperature, geometry3D, i, j, k, nDim, m, gradTx[m], gradTy[m], gradTz[m], leftOrRightIndex);
                                }
                                for (int m = 0; m < numberOfSpecies; ++ m)
                                {
                                    DXDYDZ_K(primitiveVars, geometry3D, i, j, k, nDim, m + nNSEquation, gradCx[m], gradCy[m], gradCz[m], leftOrRightIndex);
                                }

                                visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(i, j, kCellGhostLayer));
                                visTurbulence = 0.0;
                                if (wallFunctionType != WALLFUNCTION::NONE)
                                {
                                    visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(i, j, kCellGhostLayer));
                                    heatCoeff = heatTransferCoeff(i, j, k, ITT) + heatTransferCoeff(i, j, kCellGhostLayer, ITT);
                                }
                            }
                            dtdx = gradTx[ITT];
                            dtdy = gradTy[ITT];
                            dtdz = gradTz[ITT];
                            viscousCoef = visLaminar + visTurbulence;

                            //! To Compute the heat flux on wall.
                            if (nChemical == 0)
                            {
                                RDouble kSpecificHeat = (visLaminar * oPrandtlLaminar + visTurbulence * oPrandtlTurbulence) * specificHeatAtConstantPressure;
                                if (wallFunctionType == WALLFUNCTION::STANDARD)
                                {
                                    if (abs(heatCoeff) > SMALL)
                                    {
                                        kSpecificHeat += heatCoeff;
                                    }
                                }
                                RDouble qx = kSpecificHeat * dtdx;
                                RDouble qy = kSpecificHeat * dtdy;
                                RDouble qz = kSpecificHeat * dtdz;
                                //! This value is not the real heat rate on wall, it is just the first term on the right side of the energy equation.\n
                                //! To obtain the real value, hf = nxs * qx + nys * qy + nzs * qz, and hf is just the non-dimensional variables.\n
                                //! To obtain the dimensional value of heat rate, hf = (nxs * qx + nys * qy + nzs * qz) * (ref_miu * ref_vel^2/ref_len).
                                heatFlux += -(normalComponetX * qx + normalComponetY * qy + normalComponetZ * qz) * oRefReNumber;
                            }
                            else      //! Non-equilibrium gas.
                            {
                                //! Compute Heat conductivity.
                                RDouble heatConductivity[3] = { 0.0, 0.0, 0.0 };
                                //! To obtain the enthalpy of each species.
                                gas->GetEverySpeciesEnthalpy(wallTemperature[ITT], wallTemperature[ITV], wallTemperature[ITE], speciesEnthalpy);

                                //! To Modify the viscosity on wall.
                                gas->ComputeViscosityByPrimitive(wallVariables, wallTemperature[ITT], wallTemperature[ITE], nonDimensionalSutherlandTemperature, visLaminar);
                                viscousCoef = visLaminar + visTurbulence;

                                //! To obtain the mass diffusion coefficients of each species.
                                gas->GetSpeciesMassDiffusionCoef(wallVariables, visLaminar, visTurbulence, densityDiffusivity);

                                //! To compute heat conductivity of mixture gas using the Eucken formula and Wassilewa formula.
                                gas->ComputeMixtureGasHeatConductivity(wallVariables, wallTemperature[ITT], wallTemperature[ITV], wallTemperature[ITE], heatConductivity[ITT], heatConductivity[ITV], heatConductivity[ITE]);

                                if (nSurfGradMethod > 0)
                                {
                                    //! Compute the heat flux caused by the gradient of translation-rotation temperature near wall.
                                    Qtr = -heatConductivity[ITT] * (normalComponetX * gradTx[ITT] + normalComponetY * gradTy[ITT] + normalComponetZ * gradTz[ITT]) * oRefReNumber;

                                    //! Compute the heat flux caused by the gradient of species enthalpies near wall.
                                    Qs = 0.0;
                                    for (int s = 0; s < numberOfSpecies; ++ s)
                                    {
                                        //! Compute the diffusion coefficient.
                                        RDouble diffusionCoefficient = densityDiffusivity[s] * speciesEnthalpy[s];
                                        //! Compute the mass diffusion energy.
                                        Qs += -diffusionCoefficient * (normalComponetX * gradCx[s] + normalComponetY * gradCy[s] + normalComponetZ * gradCz[s]);
                                    }
                                    Qs *= oRefReNumber;

                                    //! Compute the heat flux caused by the gradient of vibration-electron temperature near wall.
                                    Qv = 0.0;
                                    Qe = 0.0;
                                    if (nTemperatureModel == 2)
                                    {
                                        Qv = -(heatConductivity[ITV] + heatConductivity[ITE]) * (normalComponetX * gradTx[ITV] + normalComponetY * gradTy[ITV] + normalComponetZ * gradTz[ITV]) * oRefReNumber;
                                    }
                                    else if (nTemperatureModel == 3)
                                    {
                                        Qv = -heatConductivity[ITV] * (normalComponetX * gradTx[ITV] + normalComponetY * gradTy[ITV] + normalComponetZ * gradTz[ITV]) * oRefReNumber;
                                        Qe = -heatConductivity[ITE] * (normalComponetX * gradTx[ITE] + normalComponetY * gradTy[ITE] + normalComponetZ * gradTz[ITE]) * oRefReNumber;
                                    }
                                }
                                else
                                {
                                    //! Compute the heat flux caused by the gradient of translation-rotation temperature near wall.
                                    Qtr = heatConductivity[ITT] * (temperature(i, j, k, ITT) - wallTemperature[ITT]) / deltaHeight * oRefReNumber;

                                    //! Compute the heat flux caused by the gradient of species enthalpies near wall.
                                    Qs = 0.0;
                                    for (int s = 0; s < numberOfSpecies; ++ s)
                                    {
                                        //! Compute the diffusion coefficient.
                                        RDouble diffusionCoefficient = densityDiffusivity[s] * speciesEnthalpy[s];
                                        //! Compute the mass diffusion energy.
                                        Qs += diffusionCoefficient * (primitiveVars(i, j, k, nNSEquation + s) - wallVariables[nNSEquation + s]) / deltaHeight;
                                    }
                                    Qs *= oRefReNumber;

                                    //! Compute the heat flux caused by the gradient of vibration-electron temperature near wall.
                                    Qv = 0.0;
                                    Qe = 0.0;
                                    if (nTemperatureModel == 2)
                                    {
                                        Qv = (heatConductivity[ITV] + heatConductivity[ITE]) * (temperature(i, j, k, ITV) - wallTemperature[ITV]) / deltaHeight * oRefReNumber;
                                    }
                                    else if (nTemperatureModel == 3)
                                    {
                                        Qv = heatConductivity[ITV] * (temperature(i, j, k, ITV) - wallTemperature[ITV]) / deltaHeight * oRefReNumber;
                                        Qe = heatConductivity[ITE] * (temperature(i, j, k, ITE) - wallTemperature[ITE]) / deltaHeight * oRefReNumber;
                                    }
                                }

                                //! Set the total heat flux on wall.
                                heatFlux = Qtr + Qs + Qv + Qe;
                                Qtr *= refHeatFluxDimension;
                                Qs *= refHeatFluxDimension;
                                Qv *= refHeatFluxDimension;
                                Qe *= refHeatFluxDimension;
                            }

                            Qw = heatFlux * refHeatFluxDimension;
                            HeatCoef = 2.0 * heatFlux;
                            RDouble deltaH = refTotalEnthalpy - wallTotalEnthalpy;
                            StantonNumber = Qw * 1000 / (refDimensionalDensity * refVelocity * deltaH);

                            RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);
                            //! Stress components.
                            RDouble txx = viscousCoef * (two * dudx - divv2p3);
                            RDouble tyy = viscousCoef * (two * dvdy - divv2p3);
                            RDouble tzz = viscousCoef * (two * dwdz - divv2p3);
                            RDouble txy = viscousCoef * (dudy + dvdx);
                            RDouble txz = viscousCoef * (dudz + dwdx);
                            RDouble tyz = viscousCoef * (dvdz + dwdy);

                            fvsx = -two * (nx * txx + ny * txy + nz * txz) * oRefReNumber;
                            fvsy = -two * (nx * txy + ny * tyy + nz * tyz) * oRefReNumber;
                            fvsz = -two * (nx * txz + ny * tyz + nz * tzz) * oRefReNumber;

                            dfx += fvsx;
                            dfy += fvsy;
                            dfz += fvsz;

                            delete geometry3D;
                        }

                        RDouble dir = cosa * fvsx + sina * fvsy;

                        RDouble cfn = fvsx * normalComponetX + fvsy * normalComponetY + fvsz * normalComponetZ;
                        RDouble cfxt = fvsx - cfn * normalComponetX;
                        RDouble cfyt = fvsy - cfn * normalComponetY;
                        RDouble cfzt = fvsz - cfn * normalComponetZ;

                        RDouble cft = sqrt(cfxt * cfxt + cfyt * cfyt + cfzt * cfzt);

                        RDouble frictionCoef = cft / surfaceArea * SIGN(1.0, dir);

                        RDouble densityWall = primitiveVars(i, j, k, IR);
                        RDouble taoWall = half * cft / surfaceArea;

                        RDouble yPlus = 0.0;
                        if (viscousType > INVISCID)
                        {
                           // RDouble3D &wallDistant = *grid->GetWallDist();
                           // RDouble wallDistantOfFirstCell = wallDistant(i, j, k);   //by dms
                            RDouble wallDistantOfFirstCell = deltaHeight;
                            yPlus = refReNumber * sqrt(taoWall * densityWall) * wallDistantOfFirstCell / visLaminar;
                        }

                        //if (nRapidFlowfield > 0)
                        //{
                        //    pressCoef = (*surfacePressure)(indexOfWall, indexOfCell, 0);
                        //    Pw = (*surfacePressure)(indexOfWall, indexOfCell, 1) * refDynamicPressure;
                        //}

                        visualVariables[VISUAL_WALL_CP] = pressCoef;
                        visualVariables[VISUAL_WALL_CF] = frictionCoef;
                        visualVariables[VISUAL_WALL_YPLUS] = yPlus;
                        visualVariables[VISUAL_WALL_QNONDIM] = heatFlux;
                        visualVariables[VISUAL_WALL_QDIM] = Qw;
                        visualVariables[VISUAL_WALL_PW] = Pw;
                        visualVariables[VISUAL_WALL_TW] = Tw;
                        visualVariables[VISUAL_WALL_RHOW] = Rhow;
                        visualVariables[VISUAL_WALL_VX] = wallVariables[IU] * refVelocity;
                        visualVariables[VISUAL_WALL_VY] = wallVariables[IV] * refVelocity;
                        visualVariables[VISUAL_WALL_VZ] = wallVariables[IW] * refVelocity;
                        visualVariables[VISUAL_WALL_VS] = absVelocity;
                        visualVariables[VISUAL_WALL_ST] = StantonNumber;
                        visualVariables[VISUAL_WALL_CH] = HeatCoef;

                        if (nChemical > 0)
                        {
                            visualVariables[VISUAL_WALL_QT] = Qtr;
                            visualVariables[VISUAL_WALL_QS] = Qs;
                            visualVariables[VISUAL_WALL_QV] = Qv;
                            visualVariables[VISUAL_WALL_QE] = Qe;
                            for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
                            {
                                visualVariables[VISUAL_SPECIES_START + iSpecies] = wallVariables[nNSEquation + iSpecies];
                            }
                        }

                        if (nSlipBCModel > 0)//The slip temperature.
                        {
                            visualVariables[VISUAL_SLIP_TS] = wallTemperature[ITT] * referenceTemperature;
                            visualVariables[VISUAL_SLIP_TV] = wallTemperature[ITV] * referenceTemperature;
                            visualVariables[VISUAL_SLIP_TE] = wallTemperature[ITE] * referenceTemperature;
                            visualVariables[VISUAL_SLIP_DTS] = deltaT * referenceTemperature;
                        }

                        int nIndex = -1;
                        for (int iVar = 0; iVar < nWallVariables; ++ iVar)
                        {
                            int varType = visualVariablesType[iVar];
                            int varIndex = variablesOrder[varType];
                            //The mass fractions of species will be exported.
                            if (varType == VISUAL_WALL_NS)
                            {
                                continue;
                            }
                            else //Obtain the other variables.
                            {
                                ++ nIndex;
                                dataStore(i, j, k, nIndex) = visualVariables[varIndex];
                            }
                        }

                        //Obtain the mass fractions of species.
                        if (nWallVariables != dataNumber)
                        {
                            for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
                            {
                                dataStore(i, j, k, nWallVariables + iSpecies - 1) = visualVariables[VISUAL_SPECIES_START + iSpecies];
                            }
                        }

                        isWall(i, j, k) = 1.0;
                        //Next grid cell.
                        ++ indexOfCell;
                    }//i end
                }//j end
            }//k end
            ++ indexOfWall; //Next surface region.
        }//iBCRegion end

        GhostCell3D(dataStore, ni, nj, nk, dataNumber);

        CommunicateInterfaceValue(grid, &isWall, "wallCell");
        CommunicateInterfaceValue(grid, &dataStore, "wallAircoefData", dataNumber);

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();

            if (!IsWall(BCType) && BCType != PHENGLEI::ABLATION_SURFACE)
            {
                continue;
            }

            int ist, ied, jst, jed, kst, ked;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            int s_nd = structBC->GetFaceDirection();
            if (s_nd == 0)
            {
                ked++ ;
                jed++ ;
            }
            else if (s_nd == 1)
            {
                ied++ ;
                ked++ ;
            }
            else
            {
                ied++ ;
                jed++ ;
            }

            if (nDim == TWO_D)
            {
                ked = 1;
            }

            int *s_lr3d = structBC->GetFaceDirectionIndex();
            int nSurface = structBC->GetFaceDirection() + 1;
            int id, jd, kd;
            GetBCFaceIDX(s_lr3d, id, jd, kd);

            vector <RDouble> vectorXCoordinate, vectorYCoordinate, vectorZCoordinate;
            vectorXCoordinate.resize(0);
            vectorYCoordinate.resize(0);
            vectorZCoordinate.resize(0);
            dataCount = 0;
            coorCount = 0;

            vector < vector <RDouble> > dumpDate;
            dumpDate.resize(dataNumber);

            kst += kd;
            ked += kd;
            jst += jd;
            jed += jd;
            ist += id;
            ied += id;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        vectorXCoordinate.push_back(xx(i, j, k));
                        vectorYCoordinate.push_back(yy(i, j, k));
                        vectorZCoordinate.push_back(zz(i, j, k));
                        coorCount++ ;

                        int iLeftCell, iRightCell, jLeftCell, jRightCell, kLeftCell, kRightCell;
                        if ((i == ist) && isWall(i - 1, j, k) < 0.5)
                        {
                            iLeftCell = i;
                        }
                        else
                        {
                            iLeftCell = i - 1;
                        }

                        if (i == ied && isWall(i, j, k) < 0.5)
                        {
                            iRightCell = i - 1;
                        }
                        else
                        {
                            iRightCell = i;
                        }

                        if (nSurface == 1)
                        {
                            iLeftCell = i - id;
                            iRightCell = i - id;
                        }

                        if (j == jst && isWall(i, j - 1, k) < 0.5)
                        {
                            jLeftCell = j;
                        }
                        else
                        {
                            jLeftCell = j - 1;
                        }

                        if (j == jed && isWall(i, j, k) < 0.5)
                        {
                            jRightCell = j - 1;
                        }
                        else
                        {
                            jRightCell = j;
                        }

                        if (nSurface == 2)
                        {
                            jLeftCell = j - jd;
                            jRightCell = j - jd;
                        }

                        if (nDim == 3)
                        {
                            if (k == kst && isWall(i, j, k - 1) < 0.5)
                            {
                                kLeftCell = k;
                            }
                            else
                            {
                                kLeftCell = k - 1;
                            }

                            if (k == ked && isWall(i, j, k) < 0.5)
                            {
                                kRightCell = k - 1;
                            }
                            else
                            {
                                kRightCell = k;
                            }

                            if (nSurface == 3)
                            {
                                kLeftCell = k - kd;
                                kRightCell = k - kd;
                            }
                        }
                        else
                        {
                            kLeftCell = 1;
                            kRightCell = 1;
                        }

                        RDouble nodeValue = 0.0;
                        for (int iData = 0; iData < dataNumber; ++ iData)
                        {
                            nodeValue = eighth * (dataStore(iLeftCell, jLeftCell, kLeftCell, iData) + dataStore(iLeftCell, jRightCell, kLeftCell, iData) +
                                dataStore(iLeftCell, jLeftCell, kRightCell, iData) + dataStore(iLeftCell, jRightCell, kRightCell, iData) +
                                dataStore(iRightCell, jLeftCell, kLeftCell, iData) + dataStore(iRightCell, jRightCell, kLeftCell, iData) +
                                dataStore(iRightCell, jLeftCell, kRightCell, iData) + dataStore(iRightCell, jRightCell, kRightCell, iData));
                            dumpDate[iData].push_back(nodeValue);
                        }
                    }
                }
            }

            TecioMission = 1;
            PHWrite(wallDistributeData, &TecioMission, 1);

            int IMxOrNumPts = ied - ist + 1;
            int JMxOrNumElements = jed - jst + 1;
            int KMxOrNumFaces = ked - kst + 1;

            PHWrite(wallDistributeData, &IMxOrNumPts, 1);
            PHWrite(wallDistributeData, &JMxOrNumElements, 1);
            PHWrite(wallDistributeData, &KMxOrNumFaces, 1);

            PHWrite(wallDistributeData, coorCount);
            PHWrite(wallDistributeData, vectorXCoordinate, coorCount);
            PHWrite(wallDistributeData, vectorYCoordinate, coorCount);
            PHWrite(wallDistributeData, vectorZCoordinate, coorCount);

            int ValueLocation = 1;
            for (int iData = 0; iData < dataNumber; ++ iData)
            {
                PHWrite(wallDistributeData, coorCount);
                PHWrite(wallDistributeData, ValueLocation);
                PHWrite(wallDistributeData, dumpDate[iData], coorCount);
            }
        }

        delete[]wallVariables;
    }

    //! To obtain the derivatives of the designate primary variable whose sequence number in the array is marked with the valued of nIndex.
    void NSSolverStruct::GetDerivatives(int nIndex, int iSurface, int ndim, int nlr, int i, int j, int k, int ni, int nj, int nk, RDouble4D &q,
        RDouble4D &xfv, RDouble4D &yfv, RDouble4D &zfv, RDouble3D &vol, RDouble &dvdx, RDouble &dvdy, RDouble &dvdz)
    {
        int p, m;
        //! Allocate.
        Geometry3D *geometry3D = new Geometry3D;

        if (iSurface == 1)
        {
            p = min(i, ni - 1);
            m = max(i - 1, 1);
            //! Initialize.
            Get_GEO_I(i, j, k, m, p, ndim, geometry3D, xfv, yfv, zfv, vol);
            //! To compute the derivatives.
            DXDYDZ_I(q, geometry3D, i, j, k, ndim, nIndex, dvdx, dvdy, dvdz, nlr);
        }
        else if (iSurface == 2)
        {
            p = min(j, nj - 1);
            m = max(j - 1, 1);
            //! Initialize.
            Get_GEO_J(i, j, k, m, p, ndim, geometry3D, xfv, yfv, zfv, vol);
            //! To compute the derivatives.
            DXDYDZ_J(q, geometry3D, i, j, k, ndim, nIndex, dvdx, dvdy, dvdz, nlr);
        }
        else
        {
            p = min(k, nk - 1);
            m = max(k - 1, 1);
            //! Initialize.
            Get_GEO_K(i, j, k, m, p, ndim, geometry3D, xfv, yfv, zfv, vol);
            //! To compute the derivatives.
            DXDYDZ_K(q, geometry3D, i, j, k, ndim, nIndex, dvdx, dvdy, dvdz, nlr);
        }

        delete geometry3D;
    }

    void NSSolverStruct::ExportSurfaceVariables(ActionKey *actkey)
    {
        int level = actkey->level;
        if (level != 0)
        {
            return;
        }
        int nDumpSurfaceInfo = GlobalDataBase::GetIntParaFromDB("nDumpSurfaceInfo");
        if (nDumpSurfaceInfo == 1)
        {
            WriteSurfaceInfo(actkey);
            return;
        }

        StructGrid *grid = StructGridCast(GetGrid(level));
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        //! The primitive temperatures.
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        //! Obtain the viscosity of turbulence.
        RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));

        Param_NSSolverStruct *parameters = GetControlParameters();

        RDouble refReNumber = parameters->GetRefReNumber();

        RDouble refGama = parameters->GetRefGama();
        RDouble oPrandtlLaminar = parameters->GetoPrandtlLaminar();
        RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

        RDouble refMachNumber = parameters->GetRefMachNumber();
        int viscousType = parameters->GetViscousType();

        int mTT = parameters->GetmTT();
        int mTV = parameters->GetmTV();
        int mTE = parameters->GetmTE();

        RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
        RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
        RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());
        //! The area normal equals to the area multiplies the normal of the face.
        RDouble4D &faceAreaNormalComponentX = *(grid->GetFaceVectorX());
        RDouble4D &faceAreaNormalComponentY = *(grid->GetFaceVectorY());
        RDouble4D &faceAreaNormalComponentZ = *(grid->GetFaceVectorZ());
        RDouble4D &faceArea = *(grid->GetFaceArea());
        RDouble3D &cellVolume = *(grid->GetCellVolume());

        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        //! Compute the number of species.
        int nSpecies = nLaminar + nChemical - nNSEquation;
        //! The total number of equation.
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

        int nSolidSurface = grid->GetNumberOfWallCell();
        if (nSolidSurface == 0)
        {
            DataContainer *cdata = actkey->GetData();
            return;
        }

        RDouble speciesEnthalpy[MAX_SPECIES_NUM] = { 0 }, densityDiffusivity[MAX_SPECIES_NUM] = { 0 };
        string *speciesName = 0;
        RDouble3D *surfaceMassFraction;
        if (nChemical > 0 && nSpecies > 0)   //! Allocate the memory.
        {
            speciesName = gas->GetNameOfSpecies();
            surfaceMassFraction = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
        }
        //! Allocate the dynamic memory space to save the primitive variables on wall.
        RDouble *wallVariables = new RDouble[nEquation]();

        //! The type of catalytic wall condition.
        RDouble catalyticCoef = parameters->GetCatalyticCoef();
        //! Temperature on wall.
        RDouble temperatureWall = parameters->GetWallTemperature();
        RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
        //! The non-dimensional temperature on wall.
        RDouble temperatureWallNonDimensional = temperatureWall / refDimensionalTemperature;
        RDouble2D *surfaceTemperature = 0;
        if (temperatureWall >= 0.0)
        {
            surfaceTemperature = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));
        }

        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        RDouble3D *surfaceSlipVariables = NULL;
        if (nSlipBCModel > 0)
        {
            surfaceSlipVariables = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceSlipVariables"));
        }

        RDouble2D *firstLayerHeight = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("firstLayerHeight"));

        //! The heating rate on wall. totalHeatViaGradients and totalHeatWall is the total heat rate, the former is computed with gradients,\n
        //! while the latter is computed with the first-order difference. transRotationHeat, vibrationHeat, eletronHeat, and speciesHeat are \n
        //! heat rate of translation-rotation, vibration, electron and species enthalpy, respectively.
        RDouble totalHeatViaGradients = 0.0, totalHeatWall = 0.0;
        RDouble transRotationHeat = 0.0, vibrationHeat = 0.0, eletronHeat = 0.0, speciesHeat = 0.0;

        //! The computations of Qhf and Qw are different. Qhf applies the partial derivative of temperature to compute heating rate,\n
        //! while Qw equals to (T1 - Tw)/dy, T1 denotes the temperature of the cell near the surface, Tw is the temperature on wall,\n
        //! and dy is computed by projecting the distant vector to the surface normal, the distant vector is constructed from the center\n
        //! of the surface to that of the cell near this surface.

        RDouble wallTemperature[3] = { temperatureWallNonDimensional, temperatureWallNonDimensional, temperatureWallNonDimensional };
        RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();

        //! Obtain the reference values.
        RDouble refLenth = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
        RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
        RDouble refVelocity = parameters->GetRefDimensionalVelocity();
        RDouble refViscosity = refDimensionalDensity * refVelocity * refLenth / refReNumber;
        RDouble refDynamicPressure = refDimensionalDensity * refVelocity * refVelocity;
        RDouble refHeatFluxDimension = refViscosity * refVelocity * refVelocity / refLenth;

        std::ostringstream oss;
        using namespace IDX;

        vector<string> title_tecplot;
        title_tecplot.push_back("title=\"Flow Fields of PHengLEI\"");
        title_tecplot.push_back("variables=");
        title_tecplot.push_back("\"x\"");
        title_tecplot.push_back("\"y\"");
        title_tecplot.push_back("\"z\"");
        title_tecplot.push_back("\"rho\""); //(kg/m3)
        title_tecplot.push_back("\"p\""); //(Pa)
        title_tecplot.push_back("\"Ttr\""); //(K)
        title_tecplot.push_back("\"Qhf\""); //(W/m2)
        title_tecplot.push_back("\"Qw\""); //(W/m2)
        title_tecplot.push_back("\"Qtr\""); //(W/m2)
        title_tecplot.push_back("\"Qv\""); //(W/m2)
        title_tecplot.push_back("\"Qe\""); //(W/m2)
        title_tecplot.push_back("\"Qs\""); //(W/m2)

        for (int iSpecies = 0; iSpecies < nSpecies; iSpecies++ )
        {
            title_tecplot.push_back(speciesName[iSpecies]);
        }

        int nDim = GetDim();
        int js = 1;
        int ks = 1;
        if (nDim == 2)
        {
            ks = 0;
        }
        if (nDim == 1)
        {
            js = 0;
        }

        //! The non-dimensional value of constant pressure specific heat.
        RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);
        //! The universal gas constant.
        RDouble gasConstant = gas->GetUniversalGasConstant();

        int nIndex = 0, indexOfWall = 0, indexOfCell = 0;
        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();    //! Obtain the number of boundary condition regions.

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int nBCType = structBC->GetBCType();

            int ist, ied, jst, jed, kst, ked;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            if (!IsWall(nBCType) && nBCType != PHENGLEI::ABLATION_SURFACE)
            {
                nIndex += (ied - ist + 1) * (jed - jst + 1) * (ked - kst + 1);
                continue;
            }

            //! Obtain the surface information.
            int *faceDirectionIndex = structBC->GetFaceDirectionIndex();

            int nSurface = structBC->GetFaceDirection() + 1;
            int leftOrRightIndex = structBC->GetFaceLeftOrRightIndex();

            RDouble normalComponetX, normalComponetY, normalComponetZ, surfaceArea, dtdx, dtdy, dtdz, dqdx, dqdy, dqdz;

            for (std::size_t i = 0; i < title_tecplot.size(); ++ i)
            {
                oss << title_tecplot[i] << "\n";
            }
            oss << "zone  i = " << ied - ist + 1
                << " j = " << jed - jst + 1
                << " k = " << ked - kst + 1
                << " f = point " << " \n";

            indexOfCell = 0;
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        //! i, j, k - is first CELL index near SOLID_SURFACE.
                        //! iWall, jWall, kWall - is NODE index on SOLID_SURFACE.
                        int iWall, jWall, kWall;
                        structBC->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);

                        //! The index of ghost cell on the negative first layer.
                        int iCellGhostLayer = i + faceDirectionIndex[0];
                        int jCellGhostLayer = j + faceDirectionIndex[1];
                        int kCellGhostLayer = k + faceDirectionIndex[2];

                        normalComponetX = faceNormalComponentX(iWall, jWall, kWall, nSurface);
                        normalComponetY = faceNormalComponentY(iWall, jWall, kWall, nSurface);
                        normalComponetZ = faceNormalComponentZ(iWall, jWall, kWall, nSurface);
                        surfaceArea = faceArea(iWall, jWall, kWall, nSurface);
                        RDouble deltaHeight = (*firstLayerHeight)(indexOfWall, indexOfCell);

                        RDouble xCenterWall, yCenterWall, zCenterWall;
                        grid->FaceCoor(iWall, jWall, kWall, nSurface, xCenterWall, yCenterWall, zCenterWall);

                        //! To obtain the real value, hf = nxs * qx + nys * qy + nzs * qz, and hf is just the non-dimensional variables.\n
                        //! To obtain the dimensional value of heat rate, hf = (nxs * qx + nys * qy + nzs * qz) * (ref_miu * ref_vel^2/ref_len).
                        //! To compute the heat rate on wall.
                        totalHeatViaGradients = 0.0;
                        totalHeatWall = 0.0;
                        transRotationHeat = 0.0;
                        vibrationHeat = 0.0;
                        eletronHeat = 0.0;
                        speciesHeat = 0.0;

                        //! Obtain the viscosity.
                        RDouble viscosityLaminar, viscosityTurbulence;
                        if (viscousType > INVISCID)
                        {
                            viscosityLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(iCellGhostLayer, jCellGhostLayer, kCellGhostLayer));
                            viscosityTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(iCellGhostLayer, jCellGhostLayer, kCellGhostLayer));
                        }
                        else
                        {
                            viscosityLaminar = 0.0;
                            viscosityTurbulence = 0.0;
                        }

                        //! To obtain the variables on wall.
                        if (nSlipBCModel > 0)
                        {
                            for (int m = 0; m < 3; ++ m)
                            {
                                wallTemperature[m] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m);
                                wallVariables[m + 1] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m + 3);
                            }
                        }
                        else
                        {
                            if (viscousType > INVISCID)
                            {
                                wallVariables[IU] = 0.0;
                                wallVariables[IV] = 0.0;
                                wallVariables[IW] = 0.0;
                            }
                            else
                            {
                                wallVariables[IU] = primitiveVariables(i, j, k, IU);
                                wallVariables[IV] = primitiveVariables(i, j, k, IV);
                                wallVariables[IW] = primitiveVariables(i, j, k, IW);
                            }

                            if (temperatureWall < 0.0)
                            {
                                wallTemperature[ITT] = temperatures(i, j, k, mTT);
                                wallTemperature[ITV] = temperatures(i, j, k, mTV);
                                wallTemperature[ITE] = temperatures(i, j, k, mTE);
                            }
                            else
                            {
                                wallTemperature[ITT] = (*surfaceTemperature)(indexOfWall, indexOfCell);
                                wallTemperature[ITV] = wallTemperature[ITT];
                                wallTemperature[ITE] = wallTemperature[ITT];
                            }
                        }

                        //! Obtain the pressure on wall.
                        wallVariables[IP] = primitiveVariables(i, j, k, IP);
                        //! Compute Heat conductivity.
                        RDouble heatConductivity[3] = { 0.0, 0.0, 0.0 };
                        //! Initialize heat conductivity of translation and rotation.
                        heatConductivity[ITT] = specificHeatAtConstantPressure * (viscosityLaminar * oPrandtlLaminar + viscosityTurbulence * oPrandtlTurbulence);
                        wallVariables[IR] = refGama * refMachNumber * refMachNumber * wallVariables[IP] / wallTemperature[ITT];

                        //! Compute the heat rate of the species diffusion, added by Li Peng on Mar 11,2019.
                        //! Chemical reaction flow.
                        if (nChemical >= 1)
                        {
                            //! The isothermal wall condition and full catalytic condition.
                            if (ABS(catalyticCoef) <= EPSILON && nBCType != PHENGLEI::ABLATION_SURFACE) //Fully non-catalytic wall.
                            {
                                for (int iLaminar = nNSEquation; iLaminar < nNSEquation + nSpecies; ++ iLaminar)
                                {
                                    wallVariables[iLaminar] = primitiveVariables(i, j, k, iLaminar);
                                }
                            }
                            else //Fully catalytic wall or finite catalytic wall.
                            {
                                for (int iLaminar = nNSEquation; iLaminar < nNSEquation + nSpecies; ++ iLaminar)
                                {
                                    wallVariables[iLaminar] = (*surfaceMassFraction)(indexOfWall, indexOfCell, iLaminar - nNSEquation);
                                }
                            }

                            // Compute density on surface.
                            RDouble massReciprocal = 0.0;
                            if (nTemperatureModel > 1)    //! Multi-Temperature model.
                            {
                                RDouble ceDivideMe = 0.0;
                                massReciprocal = gas->ComputeMolecularWeightReciprocalWithoutElectron(wallVariables, ceDivideMe);
                                wallVariables[IR] = wallVariables[IP] / (gasConstant * wallTemperature[ITT] * massReciprocal + gasConstant * wallTemperature[ITE] * ceDivideMe);
                            }
                            else     //! One-Temperature model.
                            {
                                //! Obtain molecular weight of the mixture gas.
                                massReciprocal = gas->ComputeMolecularWeightReciprocal(wallVariables);
                                wallVariables[IR] = wallVariables[IP] / (gasConstant * wallTemperature[ITT] * massReciprocal);
                            }

                            //! Compute internal energy of ghost cells.
                            if (nTemperatureModel == 2) // Two-Temperature model.
                            {
                                //! Set the internal energy of vibration-eletron on wall.
                                wallVariables[nLaminar + nChemical] = gas->GetMixedGasVibrationEnergy(wallVariables, wallTemperature[ITV]);
                                wallVariables[nLaminar + nChemical] += gas->GetMixedGasElectronEnergy(wallVariables, wallTemperature[ITV]);
                            }
                            else if (nTemperatureModel == 3)    //! Three-Temperature model.
                            {
                                //! Set the internal energy of vibration-eletron on wall.
                                wallVariables[nLaminar + nChemical] = gas->GetMixedGasVibrationEnergy(wallVariables, wallTemperature[ITV]);
                                wallVariables[nLaminar + nChemical + 1] = gas->GetMixedGasElectronEnergy(wallVariables, wallTemperature[ITE]);
                            }

                            //! To obtain the enthalpy of each species.
                            gas->GetEverySpeciesEnthalpy(wallTemperature[ITT], wallTemperature[ITV], wallTemperature[ITE], speciesEnthalpy);

                            //! To Modify the viscosity on wall.
                            gas->ComputeViscosityByPrimitive(wallVariables, wallTemperature[ITT], wallTemperature[ITE], nonDimensionalSutherlandTemperature, viscosityLaminar);

                            //! To obtain the mass diffusion coefficients of each species.
                            gas->GetSpeciesMassDiffusionCoef(wallVariables, viscosityLaminar, viscosityTurbulence, densityDiffusivity);

                            //! To compute heat conductivity of mixture gas using the Eucken formula and Wassilewa formula.
                            gas->ComputeMixtureGasHeatConductivity(wallVariables, wallTemperature[ITT], wallTemperature[ITV], wallTemperature[ITE], heatConductivity[0], heatConductivity[1], heatConductivity[2]);

                            RDouble diffusionCoefficient = 0.0;
                            dqdx = 0.0;
                            dqdy = 0.0;
                            dqdz = 0.0;
                            //! Get the partial derivatives of species.
                            for (int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies)
                            {
                                //! To compute the gradients of species density on wall.
                                RDouble dcdx, dcdy, dcdz;
                                GetDerivatives(nNSEquation + iSpecies, nSurface, nDim, leftOrRightIndex, i, j, k, ni, nj, nk, primitiveVariables,
                                    faceAreaNormalComponentX, faceAreaNormalComponentY, faceAreaNormalComponentZ, cellVolume, dcdx, dcdy, dcdz);
                                //! Compute the diffusion coefficient.
                                diffusionCoefficient = densityDiffusivity[iSpecies] * speciesEnthalpy[iSpecies];
                                dqdx += diffusionCoefficient * dcdx;
                                dqdy += diffusionCoefficient * dcdy;
                                dqdz += diffusionCoefficient * dcdz;

                                //! Compute the mass diffusion energy.
                                speciesHeat += diffusionCoefficient * (primitiveVariables(i, j, k, nNSEquation + iSpecies) - wallVariables[nNSEquation + iSpecies]) / deltaHeight;
                            }

                            //! The non-dimensional value of heat rate, to obtain the dimensional value by multiplying the Reynold number.
                            totalHeatViaGradients += -leftOrRightIndex * (normalComponetX * dqdx + normalComponetY * dqdy + normalComponetZ * dqdz); //The species diffusion energy

                            //! The first method using the difference of the mass fractions in the actual cell and the ghost cell.
                            totalHeatWall += speciesHeat;
                        }

                        //! Compute the heat rate caused by the gradient of temperature.
                        if (viscousType > INVISCID)
                        {
                            //! Obtain the partial derivative of temperature.
                            GetDerivatives(ITT, nSurface, nDim, leftOrRightIndex, i, j, k, ni, nj, nk, temperatures,
                                faceAreaNormalComponentX, faceAreaNormalComponentY, faceAreaNormalComponentZ, cellVolume, dtdx, dtdy, dtdz);
                            dqdx = heatConductivity[0] * dtdx;
                            dqdy = heatConductivity[0] * dtdy;
                            dqdz = heatConductivity[0] * dtdz;
                            //! The non-dimensional value of heat rate, to obtain the dimensional value by multiplying the Reynold number.
                            totalHeatViaGradients += -leftOrRightIndex * (normalComponetX * dqdx + normalComponetY * dqdy + normalComponetZ * dqdz);     //! Used for comparison.

                            //! Compute heat rate on wall caused by gradient of temperature.
                            transRotationHeat = heatConductivity[0] * (temperatures(i, j, k, ITT) - wallTemperature[ITT]) / deltaHeight;
                            totalHeatWall += transRotationHeat;

                            if (nTemperatureModel == 2)     //! Two-Temperature model.
                            {
                                GetDerivatives(ITV, nSurface, nDim, leftOrRightIndex, i, j, k, ni, nj, nk, temperatures,
                                    faceAreaNormalComponentX, faceAreaNormalComponentY, faceAreaNormalComponentZ, cellVolume, dtdx, dtdy, dtdz);
                                dqdx = (heatConductivity[ITV] + heatConductivity[ITE]) * dtdx;
                                dqdy = (heatConductivity[ITV] + heatConductivity[ITE]) * dtdy;
                                dqdz = (heatConductivity[ITV] + heatConductivity[ITE]) * dtdz;
                                totalHeatViaGradients += -leftOrRightIndex * (normalComponetX * dqdx + normalComponetY * dqdy + normalComponetZ * dqdz); //! Addition of vibration-electron energy.

                                vibrationHeat = (heatConductivity[ITV] + heatConductivity[ITE]) * (temperatures(i, j, k, ITV) - wallTemperature[ITV]) / deltaHeight;
                                totalHeatWall += vibrationHeat;
                            }
                            else if (nTemperatureModel == 3)    //! Three-Temperature model.
                            {
                                GetDerivatives(ITV, nSurface, nDim, leftOrRightIndex, i, j, k, ni, nj, nk, temperatures,
                                    faceAreaNormalComponentX, faceAreaNormalComponentY, faceAreaNormalComponentZ, cellVolume, dtdx, dtdy, dtdz);
                                dqdx = heatConductivity[ITV] * dtdx;
                                dqdy = heatConductivity[ITV] * dtdy;
                                dqdz = heatConductivity[ITV] * dtdz;
                                totalHeatViaGradients += -leftOrRightIndex * (normalComponetX * dqdx + normalComponetY * dqdy + normalComponetZ * dqdz);// Addition of vibration energy.

                                GetDerivatives(ITE, nSurface, nDim, leftOrRightIndex, i, j, k, ni, nj, nk, temperatures,
                                    faceAreaNormalComponentX, faceAreaNormalComponentY, faceAreaNormalComponentZ, cellVolume, dtdx, dtdy, dtdz);
                                dqdx = heatConductivity[ITE] * dtdx;
                                dqdy = heatConductivity[ITE] * dtdy;
                                dqdz = heatConductivity[ITE] * dtdz;
                                totalHeatViaGradients += -leftOrRightIndex * (normalComponetX * dqdx + normalComponetY * dqdy + normalComponetZ * dqdz); // Addition of electron energy.

                                vibrationHeat = heatConductivity[ITV] * (temperatures(i, j, k, ITV) - wallTemperature[ITV]) / deltaHeight;
                                eletronHeat = heatConductivity[ITE] * (temperatures(i, j, k, ITE) - wallTemperature[ITE]) / deltaHeight;
                                totalHeatWall += vibrationHeat + eletronHeat;
                            }
                        }

                        //! Trans the value of the heat rate to the dimensional value.
                        totalHeatViaGradients *= refHeatFluxDimension;
                        totalHeatWall *= refHeatFluxDimension;
                        transRotationHeat *= refHeatFluxDimension;
                        vibrationHeat *= refHeatFluxDimension;
                        eletronHeat *= refHeatFluxDimension;
                        speciesHeat *= refHeatFluxDimension;
                        //! Trans the non-dimensional values to the dimensional values.
                        wallVariables[IR] *= refDimensionalDensity;
                        wallVariables[IP] *= refDynamicPressure;
                        wallTemperature[ITT] *= refDimensionalTemperature;

                        //! Write to file.
                        int wordwidth = 20;
                        oss << setiosflags(ios::left);
                        oss << setiosflags(ios::scientific);
                        oss << setprecision(10);
                        oss << setw(wordwidth) << xCenterWall;
                        oss << setw(wordwidth) << yCenterWall;
                        oss << setw(wordwidth) << zCenterWall;
                        oss << setw(wordwidth) << wallVariables[IR];
                        oss << setw(wordwidth) << wallVariables[IP];
                        oss << setw(wordwidth) << wallTemperature[ITT];
                        oss << setw(wordwidth) << totalHeatViaGradients;
                        oss << setw(wordwidth) << totalHeatWall;
                        oss << setw(wordwidth) << transRotationHeat;
                        oss << setw(wordwidth) << vibrationHeat;
                        oss << setw(wordwidth) << eletronHeat;
                        oss << setw(wordwidth) << speciesHeat;

                        for (int iSpecies = 0; iSpecies < nSpecies; ++ iSpecies)
                        {
                            oss << setw(wordwidth) << wallVariables[iSpecies + nNSEquation];
                        }
                        oss << "\n";
                        ++ nIndex;
                        ++ indexOfCell; //Next grid cell.
                    }
                }
            }
            ++ indexOfWall; //Next surface region.
        }

        //! Deallocate the dynamic memory space.
        delete[]wallVariables;

        string str = oss.str();

        DataContainer *chemicalSurfaceData = actkey->GetData();
        chemicalSurfaceData->moveto_begin();
        chemicalSurfaceData->write(const_cast <char *> (str.data()), str.size() * sizeof(char));
    }

    void NSSolverStruct::GetSurfaceHeatingChange(ActionKey *actkey, RDouble &localMaxHeatChange)
    {
        using namespace IDX;
        int nDim = GetDim();
        StructGrid *grid = StructGridCast(GetGrid(actkey->level));

        RDouble4D &primitiveVars = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperature = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));

        Param_NSSolverStruct *parameters = GetControlParameters();
        RDouble oRefReNumber = parameters->GetoRefReNumber();
        RDouble oPrandtlLaminar = parameters->GetoPrandtlLaminar();
        RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();
        RDouble refGama = parameters->GetRefGama();
        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);

        int viscousType = parameters->GetViscousType();
        int mTT = parameters->GetmTT();
        int mTV = parameters->GetmTV();
        int mTE = parameters->GetmTE();

        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());

        RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
        RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
        RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());
        RDouble4D &faceArea = *(grid->GetFaceArea());
        RDouble3D &cellVolume = *(grid->GetCellVolume());

        RDouble gasConstant = gas->GetUniversalGasConstant();//! The universal gas constant.
        RDouble catalyticCoef = parameters->GetCatalyticCoef();//! The type of catalytic wall condition.
        RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();
        RDouble Twall = parameters->GetWallTemperature();    //! Temperature on wall.
        RDouble referenceTemperature = parameters->GetRefDimensionalTemperature();
        int nChemical = parameters->GetChemicalFlag();
        int nNSEquation = parameters->GetNSEquationNumber();
        int numberOfSpecies = parameters->GetNumberOfSpecies();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nNSEquation + numberOfSpecies + nTemperatureModel - 1;
        int wallFunctionType = GlobalDataBase::GetIntParaFromDB("wallFunctionType");

        RDouble2D *surfaceTemperature = NULL;
        if (Twall >= 0.0)
        {
            surfaceTemperature = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));
        }

        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        RDouble3D *surfaceSlipVariables = NULL;
        if (nSlipBCModel > 0)
        {
            surfaceSlipVariables = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceSlipVariables"));
        }

        RDouble2D *firstLayerHeight = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("firstLayerHeight"));
        RDouble4D &heatTransferCoeff = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("heatTransferCoeff"));

        //! Allocate the dynamic memory space to save the primitive variables on wall.
        RDouble *wallVariables = new RDouble[nEquation];
        RDouble nondimTwall = Twall / referenceTemperature;
        RDouble wallTemperature[3] = { nondimTwall, nondimTwall, nondimTwall };
        RDouble Qw = 0.0, Qtr = 0.0, Qs = 0.0, Qv = 0.0, Qe = 0.0, prevHeat = 0.0, deltaHeat = 0.0;
        RDouble maxHeatChangeRatio = 0.0, localHeatChangeRatio = 0.0;

        int nSurfGradMethod = 0;
        if (GlobalDataBase::IsExist("nSurfGradMethod", PHINT, 1))
        {
            nSurfGradMethod = GlobalDataBase::GetIntParaFromDB("nSurfGradMethod");
        }
        RDouble gradTx[3] = {0.0, 0.0, 0.0}, gradTy[3] = {0.0, 0.0, 0.0}, gradTz[3] = {0.0, 0.0, 0.0};
        RDouble gradCx[MAX_SPECIES_NUM] = {0.0}, gradCy[MAX_SPECIES_NUM] = {0.0}, gradCz[MAX_SPECIES_NUM] = {0.0};

        //! The mass fraction of oxygen and nitrogen under the fully catalytic condition.
        RDouble speciesEnthalpy[MAX_SPECIES_NUM], densityDiffusivity[MAX_SPECIES_NUM];
        string *speciesName = NULL;
        RDouble3D *surfaceMassFraction = NULL;
        if (nChemical > 0 && numberOfSpecies > 0)
        {
            speciesName = gas->GetNameOfSpecies();
            surfaceMassFraction = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
        }

        RDouble3D *surfaceHeatFlux = NULL;
        int nSurfHeatMonitor = parameters->GetSurfaceHeatingMonitor();
        if (nSurfHeatMonitor > 0)
        {
            surfaceHeatFlux = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceHeatFlux"));
        }

        //! To obtain the reference values.
        RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
        RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
        RDouble refVelocity = refDimensionalSonicSpeed * refMachNumber;
        RDouble refSquareVelocity = refVelocity * refVelocity;
        RDouble refDynamicPressure = refDimensionalDensity * refSquareVelocity;
        RDouble refEnergy = refDimensionalDensity * refVelocity * refSquareVelocity;
        RDouble refHeatFluxDimension = refEnergy * 0.001; // The dimensional value of total heat flux (kW/m2).

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        //! Obtain the number of boundary condition regions.
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        int indexOfWall = 0, indexOfCell = 0;

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();

            if (!IsWall(BCType) && BCType != PHENGLEI::ABLATION_SURFACE)
            {
                continue;
            }

            int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            //! Obtain the surface information.
            int *s_lr3d = structBC->GetFaceDirectionIndex();
            int nSurface = structBC->GetFaceDirection() + 1;
            int leftOrRightIndex = structBC->GetFaceLeftOrRightIndex();

            indexOfCell = 0;
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        //! i, j, k - is first CELL index near SOLID_SURFACE.
                        //! iWall, jWall, kWall - is NODE index on SOLID_SURFACE.
                        int iWall, jWall, kWall;
                        structBC->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);

                        RDouble normalComponetX = leftOrRightIndex * faceNormalComponentX(iWall, jWall, kWall, nSurface);
                        RDouble normalComponetY = leftOrRightIndex * faceNormalComponentY(iWall, jWall, kWall, nSurface);
                        RDouble normalComponetZ = leftOrRightIndex * faceNormalComponentZ(iWall, jWall, kWall, nSurface);
                        RDouble surfaceArea = faceArea(iWall, jWall, kWall, nSurface);

                        //! The projection of the dy called the distance between the center of the cell and that of the surface on normal vector.
                        //RDouble deltaHeight = fabs((xCellCenter - xCenterWall) * normalComponetX + (yCellCenter - yCenterWall) * normalComponetY + (zCellCenter - zCenterWall) * normalComponetZ);
                        RDouble deltaHeight = (*firstLayerHeight)(indexOfWall, indexOfCell);

                        //! The index of  first ghost cell layer near SOLID_SURFACE.
                        int iCellGhostLayer = i + s_lr3d[0];
                        int jCellGhostLayer = j + s_lr3d[1];
                        int kCellGhostLayer = k + s_lr3d[2];

                        //! To obtain the variables on wall.
                        wallVariables[IP] = primitiveVars(i, j, k, IP);
                        if (nSlipBCModel > 0)
                        {
                            for (int m = 0; m < 3; ++ m)
                            {
                                wallTemperature[m] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m);
                                wallVariables[m + 1] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m + 3);
                            }
                        }
                        else
                        {
                            if (viscousType > INVISCID)
                            {
                                wallVariables[IU] = 0.0;
                                wallVariables[IV] = 0.0;
                                wallVariables[IW] = 0.0;
                            }
                            else
                            {
                                //The velocity is used to extract the wall streamline, so this velocity is not the real velocity on wall.
                                wallVariables[IU] = primitiveVars(i, j, k, IU);
                                wallVariables[IV] = primitiveVars(i, j, k, IV);
                                wallVariables[IW] = primitiveVars(i, j, k, IW);
                            }

                            //! Set temperature on wall.
                            if (Twall < 0.0) //Adabatic wall.
                            {
                                wallTemperature[ITT] = temperature(i, j, k, mTT);
                                wallTemperature[ITV] = temperature(i, j, k, mTV);
                                wallTemperature[ITE] = temperature(i, j, k, mTE);
                            }
                            else //Isothermal wall or radiation equilibrium temperature wall.
                            {
                                wallTemperature[ITT] = (*surfaceTemperature)(indexOfWall, indexOfCell);
                                wallTemperature[ITV] = wallTemperature[ITT];
                                wallTemperature[ITE] = wallTemperature[ITT];
                            }
                        }

                        //! Compute the density on wall using the state equation.
                        wallVariables[IR] = refGama * refMachNumber * refMachNumber * wallVariables[IP] / wallTemperature[ITT];
                        if (nChemical >= 1)
                        {
                            //! To obtain the mass fraction on wall.
                            if (ABS(catalyticCoef) <= EPSILON && BCType != PHENGLEI::ABLATION_SURFACE) // Fully non-catalytic wall.
                            {
                                for (int s = nNSEquation; s < numberOfSpecies + nNSEquation; ++ s)
                                {
                                    wallVariables[s] = primitiveVars(i, j, k, s);
                                }
                            }
                            else // Fully catalytic wall or finite catalytic wall.
                            {
                                for (int s = nNSEquation; s < numberOfSpecies + nNSEquation; ++ s)
                                {
                                    wallVariables[s] = (*surfaceMassFraction)(indexOfWall, indexOfCell, s - nNSEquation);
                                }
                            }

                            RDouble massReciprocal = 0.0;
                            if (nTemperatureModel > 1)     //! Multi-Temperature model.
                            {
                                RDouble ceDivideMe = 0.0;
                                massReciprocal = gas->ComputeMolecularWeightReciprocalWithoutElectron(wallVariables, ceDivideMe);
                                wallVariables[IR] = wallVariables[IP] / (gasConstant * wallTemperature[ITT] * massReciprocal + gasConstant * wallTemperature[ITE] * ceDivideMe);
                            }
                            else         //! One-Temperature model.
                            {
                                //! Obtain molecular weight of the mixture gas.
                                massReciprocal = gas->ComputeMolecularWeightReciprocal(wallVariables);
                                wallVariables[IR] = wallVariables[IP] / (gasConstant * wallTemperature[ITT] * massReciprocal);
                            }
                        }

                        RDouble dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
                        RDouble dtdx, dtdy, dtdz;
                        RDouble visLaminar = 0.0, visTurbulence = 0.0;
                        RDouble heatCoeff = 0.0, heatFlux = 0.0;
                        if (viscousType > INVISCID)
                        {
                            Geometry3D *geometry3D = new Geometry3D;

                            if (nSurface == 1)
                            {
                                //! Get i-th face vector on the wall surface.
                                Get_GEO_I(i, j, k, iWall, iWall, nDim, geometry3D, xfv, yfv, zfv, cellVolume);
                                //! Get derivative to the wall in i-direction of u,v,w,T.
                                DXDYDZ_I(primitiveVars, geometry3D, i, j, k, nDim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                DXDYDZ_I(primitiveVars, geometry3D, i, j, k, nDim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                DXDYDZ_I(primitiveVars, geometry3D, i, j, k, nDim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                //DXDYDZ_I(temperature, geometry3D, i, j, k, nDim, ITT, dtdx, dtdy, dtdz, leftOrRightIndex);
                                for (int m = 0; m < nTemperatureModel; ++ m)
                                {
                                    DXDYDZ_I(temperature, geometry3D, i, j, k, nDim, m, gradTx[m], gradTy[m], gradTz[m], leftOrRightIndex);
                                }
                                for (int m = 0; m < numberOfSpecies; ++ m)
                                {
                                    DXDYDZ_I(primitiveVars, geometry3D, i, j, k, nDim, m + nNSEquation, gradCx[m], gradCy[m], gradCz[m], leftOrRightIndex);
                                }

                                visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(iCellGhostLayer, j, k));
                                visTurbulence = 0.0;
                                if (wallFunctionType != WALLFUNCTION::NONE)
                                {
                                    visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(iCellGhostLayer, j, k));
                                    heatCoeff = heatTransferCoeff(i, j, k, ITT) + heatTransferCoeff(iCellGhostLayer, j, k, ITT);
                                }
                            }
                            else if (nSurface == 2)
                            {
                                //! Get j-th face vector on the wall surface;.
                                Get_GEO_J(i, j, k, jWall, jWall, nDim, geometry3D, xfv, yfv, zfv, cellVolume);
                                //! Get derivative to the wall in j-direction of u,v,w,T.
                                DXDYDZ_J(primitiveVars, geometry3D, i, j, k, nDim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                DXDYDZ_J(primitiveVars, geometry3D, i, j, k, nDim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                DXDYDZ_J(primitiveVars, geometry3D, i, j, k, nDim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                //DXDYDZ_J(temperature, geometry3D, i, j, k, nDim, ITT, dtdx, dtdy, dtdz, leftOrRightIndex);
                                for (int m = 0; m < nTemperatureModel; ++ m)
                                {
                                    DXDYDZ_J(temperature, geometry3D, i, j, k, nDim, m, gradTx[m], gradTy[m], gradTz[m], leftOrRightIndex);
                                }
                                for (int m = 0; m < numberOfSpecies; ++ m)
                                {
                                    DXDYDZ_J(primitiveVars, geometry3D, i, j, k, nDim, m + nNSEquation, gradCx[m], gradCy[m], gradCz[m], leftOrRightIndex);
                                }

                                visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(i, jCellGhostLayer, k));
                                visTurbulence = 0.0;
                                if (wallFunctionType != WALLFUNCTION::NONE)
                                {
                                    visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(i, jCellGhostLayer, k));
                                    heatCoeff = heatTransferCoeff(i, j, k, ITT) + heatTransferCoeff(i, jCellGhostLayer, k, ITT);
                                }
                            }
                            else
                            {
                                //! Get k-th face vector on the wall surface.
                                Get_GEO_K(i, j, k, kWall, kWall, nDim, geometry3D, xfv, yfv, zfv, cellVolume);
                                //! Get derivative to the wall in k-direction of u,v,w,T.
                                DXDYDZ_K(primitiveVars, geometry3D, i, j, k, nDim, IU, dudx, dudy, dudz, leftOrRightIndex);
                                DXDYDZ_K(primitiveVars, geometry3D, i, j, k, nDim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                                DXDYDZ_K(primitiveVars, geometry3D, i, j, k, nDim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                                //DXDYDZ_K(temperature, geometry3D, i, j, k, nDim, ITT, dtdx, dtdy, dtdz, leftOrRightIndex);
                                for (int m = 0; m < nTemperatureModel; ++ m)
                                {
                                    DXDYDZ_K(temperature, geometry3D, i, j, k, nDim, m, gradTx[m], gradTy[m], gradTz[m], leftOrRightIndex);
                                }
                                for (int m = 0; m < numberOfSpecies; ++ m)
                                {
                                    DXDYDZ_K(primitiveVars, geometry3D, i, j, k, nDim, m + nNSEquation, gradCx[m], gradCy[m], gradCz[m], leftOrRightIndex);
                                }

                                visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(i, j, kCellGhostLayer));
                                visTurbulence = 0.0;
                                if (wallFunctionType != WALLFUNCTION::NONE)
                                {
                                    visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(i, j, kCellGhostLayer));
                                    heatCoeff = heatTransferCoeff(i, j, k, ITT) + heatTransferCoeff(i, j, kCellGhostLayer, ITT);
                                }
                            }
                            dtdx = gradTx[ITT];
                            dtdy = gradTy[ITT];
                            dtdz = gradTz[ITT];
                            RDouble viscousCoef = visLaminar + visTurbulence;

                            //! To Compute the heat flux on wall.
                            if (nChemical == 0)
                            {
                                RDouble kSpecificHeat = (visLaminar * oPrandtlLaminar + visTurbulence * oPrandtlTurbulence) * specificHeatAtConstantPressure;
                                if (wallFunctionType == WALLFUNCTION::STANDARD)
                                {
                                    if (abs(heatCoeff) > SMALL)
                                    {
                                        kSpecificHeat = heatCoeff;
                                    }
                                }
                                RDouble qx = kSpecificHeat * dtdx;
                                RDouble qy = kSpecificHeat * dtdy;
                                RDouble qz = kSpecificHeat * dtdz;
                                //! This value is not the real heat rate on wall, it is just the first term on the right side of the energy equation.\n
                                //! To obtain the real value, hf = nxs * qx + nys * qy + nzs * qz, and hf is just the non-dimensional variables.\n
                                //! To obtain the dimensional value of heat rate, hf = (nxs * qx + nys * qy + nzs * qz) * (ref_miu * ref_vel^2/ref_len).
                                heatFlux += -(normalComponetX * qx + normalComponetY * qy + normalComponetZ * qz) * oRefReNumber;
                            }
                            else      //! Non-equilibrium gas.
                            {
                                //! Compute Heat conductivity.
                                RDouble heatConductivity[3] = { 0.0, 0.0, 0.0 };
                                //! To obtain the enthalpy of each species.
                                gas->GetEverySpeciesEnthalpy(wallTemperature[ITT], wallTemperature[ITV], wallTemperature[ITE], speciesEnthalpy);

                                //! To Modify the viscosity on wall.
                                gas->ComputeViscosityByPrimitive(wallVariables, wallTemperature[ITT], wallTemperature[ITE], nonDimensionalSutherlandTemperature, visLaminar);
                                viscousCoef = visLaminar + visTurbulence;

                                //! To obtain the mass diffusion coefficients of each species.
                                gas->GetSpeciesMassDiffusionCoef(wallVariables, visLaminar, visTurbulence, densityDiffusivity);

                                //! To compute heat conductivity of mixture gas using the Eucken formula and Wassilewa formula.
                                gas->ComputeMixtureGasHeatConductivity(wallVariables, wallTemperature[ITT], wallTemperature[ITV], wallTemperature[ITE], heatConductivity[ITT], heatConductivity[ITV], heatConductivity[ITE]);

                                if (nSurfGradMethod > 0)
                                {
                                    //! Compute the heat flux caused by the gradient of translation-rotation temperature near wall.
                                    Qtr = -heatConductivity[ITT] * (normalComponetX * gradTx[ITT] + normalComponetY * gradTy[ITT] + normalComponetZ * gradTz[ITT]) * oRefReNumber;

                                    //! Compute the heat flux caused by the gradient of species enthalpies near wall.
                                    Qs = 0.0;
                                    for (int s = 0; s < numberOfSpecies; ++ s)
                                    {
                                        //! Compute the diffusion coefficient.
                                        RDouble diffusionCoefficient = densityDiffusivity[s] * speciesEnthalpy[s];
                                        //! Compute the mass diffusion energy.
                                        Qs += -diffusionCoefficient * (normalComponetX * gradCx[s] + normalComponetY * gradCy[s] + normalComponetZ * gradCz[s]);
                                    }
                                    Qs *= oRefReNumber;

                                    //! Compute the heat flux caused by the gradient of vibration-electron temperature near wall.
                                    Qv = 0.0;
                                    Qe = 0.0;
                                    if (nTemperatureModel == 2)
                                    {
                                        Qv = -(heatConductivity[ITV] + heatConductivity[ITE]) * (normalComponetX * gradTx[ITV] + normalComponetY * gradTy[ITV] + normalComponetZ * gradTz[ITV]) * oRefReNumber;
                                    }
                                    else if (nTemperatureModel == 3)
                                    {
                                        Qv = -heatConductivity[ITV] * (normalComponetX * gradTx[ITV] + normalComponetY * gradTy[ITV] + normalComponetZ * gradTz[ITV]) * oRefReNumber;
                                        Qe = -heatConductivity[ITE] * (normalComponetX * gradTx[ITE] + normalComponetY * gradTy[ITE] + normalComponetZ * gradTz[ITE]) * oRefReNumber;
                                    }
                                }
                                else
                                {
                                    //! Compute the heat flux caused by the gradient of translation-rotation temperature near wall.
                                    Qtr = heatConductivity[ITT] * (temperature(i, j, k, ITT) - wallTemperature[ITT]) / deltaHeight * oRefReNumber;

                                    //! Compute the heat flux caused by the gradient of species enthalpies near wall.
                                    Qs = 0.0;
                                    for (int s = 0; s < numberOfSpecies; ++ s)
                                    {
                                        //! Compute the diffusion coefficient.
                                        RDouble diffusionCoefficient = densityDiffusivity[s] * speciesEnthalpy[s];
                                        //! Compute the mass diffusion energy.
                                        Qs += diffusionCoefficient * (primitiveVars(i, j, k, nNSEquation + s) - wallVariables[nNSEquation + s]) / deltaHeight;
                                    }
                                    Qs *= oRefReNumber;

                                    //! Compute the heat flux caused by the gradient of vibration-electron temperature near wall.
                                    Qv = 0.0;
                                    Qe = 0.0;
                                    if (nTemperatureModel == 2)
                                    {
                                        Qv = (heatConductivity[ITV] + heatConductivity[ITE]) * (temperature(i, j, k, ITV) - wallTemperature[ITV]) / deltaHeight * oRefReNumber;
                                    }
                                    else if (nTemperatureModel == 3)
                                    {
                                        Qv = heatConductivity[ITV] * (temperature(i, j, k, ITV) - wallTemperature[ITV]) / deltaHeight * oRefReNumber;
                                        Qe = heatConductivity[ITE] * (temperature(i, j, k, ITE) - wallTemperature[ITE]) / deltaHeight * oRefReNumber;
                                    }
                                }

                                //! Set the total heat flux on wall.
                                heatFlux = Qtr + Qs + Qv + Qe;
                                Qtr *= refHeatFluxDimension;
                                Qs *= refHeatFluxDimension;
                                Qv *= refHeatFluxDimension;
                                Qe *= refHeatFluxDimension;
                            }
                            //The dimensional value of total surface heating ratio.
                            Qw = heatFlux * refHeatFluxDimension;
                        }

                        prevHeat = (*surfaceHeatFlux)(indexOfWall, indexOfCell, 0);
                        deltaHeight = heatFlux - prevHeat;
                        if (fabs(prevHeat) <= EPSILON)
                        {
                            localHeatChangeRatio = fabs(deltaHeight / EPSILON);
                        }
                        else
                        {
                            localHeatChangeRatio = fabs(deltaHeight / prevHeat);
                        }
                        (*surfaceHeatFlux)(indexOfWall, indexOfCell, 0) = heatFlux;
                        (*surfaceHeatFlux)(indexOfWall, indexOfCell, 1) = deltaHeight;

                        maxHeatChangeRatio = MAX(maxHeatChangeRatio, localHeatChangeRatio);

                        //Next grid cell.
                        ++ indexOfCell;
                    }//i end
                }//j end
            }//k end
            ++ indexOfWall; //Next surface region.
        }//iBCRegion end

        localMaxHeatChange = maxHeatChangeRatio;
    }

    void NSSolverStruct::WriteSurfaceInfo(ActionKey *actkey)
    {
        StructGrid *grid = StructGridCast(GetGrid(actkey->level));
        int nSurfaceCell = grid->GetNumberOfWallCell();
        if (nSurfaceCell == 0)
        {
            DataContainer *cdata = actkey->GetData();
            return;
        }

        RDouble4D &faceArea = *(grid->GetFaceArea());

        Param_NSSolverStruct *parameters = GetControlParameters();
        RDouble temperatureWall = parameters->GetWallTemperature();
        RDouble2D *surfaceTemperature = NULL;
        if (temperatureWall >= 0.0)
        {
            surfaceTemperature = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));
        }
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble refTemperature = parameters->GetRefDimensionalTemperature();
        RDouble refLenth = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
        RDouble refAera = refLenth * refLenth;

        std::ostringstream oss;
        using namespace IDX;
        vector<string> title_tecplot;
        title_tecplot.push_back("title=\"Surface Information of PHengLEI\"");
        title_tecplot.push_back("variables=");
        title_tecplot.push_back("\"x\"");
        title_tecplot.push_back("\"y\"");
        title_tecplot.push_back("\"z\"");
        title_tecplot.push_back("\"T\"");
        title_tecplot.push_back("\"S\"");

        int nDim = GetDim();
        int js = 1;
        int ks = 1;
        if (nDim == 2)
        {
            ks = 0;
            refAera = refLenth;
        }
        if (nDim == 1)
        {
            js = 0;
        }

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        //! Obtain the number of boundary condition regions.
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        int nIndex = 0, indexOfWall = 0, indexOfCell = 0;
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int nBCType = structBC->GetBCType();

            int ist, ied, jst, jed, kst, ked;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            if (!IsWall(nBCType))
            {
                nIndex += (ied - ist + 1) * (jed - jst + 1) * (ked - kst + 1);
                continue;
            }

            //! Obtain the surface information.
            int nSurface = structBC->GetFaceDirection() + 1;

            //! Write the title information.
            for (std::size_t i = 0; i < title_tecplot.size(); ++ i)
            {
                oss << title_tecplot[i] << "\n";
            }
            oss << "zone  i = " << ied - ist + 1
                << " j = " << jed - jst + 1
                << " k = " << ked - kst + 1
                << " f = point " << " \n";

            RDouble Tw = 1.0, surfaceAera = SMALL;
            indexOfCell = 0;
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        //! i, j, k - is first CELL index near SOLID_SURFACE.
                        //! iWall, jWall, kWall - is NODE index on SOLID_SURFACE.
                        int iWall, jWall, kWall;
                        structBC->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);

                        //Obtain the surface aera.
                        surfaceAera = faceArea(iWall, jWall, kWall, nSurface);
                        RDouble xCenterWall, yCenterWall, zCenterWall;
                        grid->FaceCoor(iWall, jWall, kWall, nSurface, xCenterWall, yCenterWall, zCenterWall);

                        //Obtain the temperature on wall.
                        if (temperatureWall >= 0.0)
                        {
                            Tw = (*surfaceTemperature)(indexOfWall, indexOfCell);
                        }
                        else
                        {
                            Tw = temperatures(i, j, k, ITT);
                        }
                        Tw *= refTemperature;
                        surfaceAera *= refAera;

                        //! Write to file.
                        int wordwidth = 20;
                        oss << setiosflags(ios::left);
                        oss << setiosflags(ios::scientific);
                        oss << setprecision(10);
                        oss << setw(wordwidth) << xCenterWall;
                        oss << setw(wordwidth) << yCenterWall;
                        oss << setw(wordwidth) << zCenterWall;
                        oss << setw(wordwidth) << Tw;
                        oss << setw(wordwidth) << surfaceAera;
                        oss << "\n";

                        ++ nIndex;
                        ++ indexOfCell; //Next grid cell.
                    }
                }
            }
            ++ indexOfWall; //Next surface region.
        }

        string str = oss.str();
        DataContainer *chemicalSurfaceData = actkey->GetData();
        chemicalSurfaceData->moveto_begin();
        chemicalSurfaceData->write(const_cast <char *> (str.data()), str.size() * sizeof(char));
    }

    void NSSolverStruct::GetResidual(ActionKey *actkey)
    {
        int level = actkey->level;
        int nEquation = GetNumberOfEquations();

        StructGrid *grid = StructGridCast(GetGrid(level));

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));

        PHSPACE::ComputeResidualonGrid(grid, actkey, residual, nEquation);
    }

    void NSSolverStruct::GetResidual(ActionKey *actkey, RDouble &localMaxRes)
    {
        int level = actkey->level;
        Param_NSSolverStruct *parameters = GetControlParameters();
        StructGrid *grid = StructGridCast(GetGrid(level));

        int varIndex = 0;
        if (GlobalDataBase::IsExist("nKeyVariableIndex", PHINT, 1))
        {
            varIndex = GlobalDataBase::GetIntParaFromDB("nKeyVariableIndex");
        }

        RDouble *maxResidualVariation = reinterpret_cast <RDouble *> (grid->GetDataPtr("maxResidualVariation"));

        localMaxRes = maxResidualVariation[varIndex];
    }

    void NSSolverStruct::GetResidual(Grid *gridIn, RDouble &localMaxRes)
    {
        Param_NSSolverStruct *parameters = GetControlParameters();
        StructGrid *grid = StructGridCast(gridIn);
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        //if (nChemical == 0)
        //{
        //    varIndex = 4; //The pressure is used.
        //}
        //else
        //{
        //    varIndex = 1; //The translation temperature is used.
        //    if (nTemperatureModel >= 2)
        //    {
        //        varIndex = 2; //The vibration temperature is used.
        //    }
        //}
        int varIndex = 0;
        if (GlobalDataBase::IsExist("nKeyVariableIndex", PHINT, 1))
        {
            varIndex = GlobalDataBase::GetIntParaFromDB("nKeyVariableIndex");
        }

        RDouble *maxResidualVariation = reinterpret_cast <RDouble *> (grid->GetDataPtr("maxResidualVariation"));

        localMaxRes = maxResidualVariation[varIndex];
    }

    void NSSolverStruct::ComputeVorticitybyQCriteria(Grid *gridIn, RDouble3D *&XVorticity, RDouble3D *&YVorticity, RDouble3D *&ZVorticity, RDouble3D *&vorticityMagnitude, RDouble3D *&strainRate, RDouble3D *&QCriteria)
    {
        StructGrid *grid = StructGridCast(gridIn);

        using namespace IDX;

        RDouble4D &gradUVWTCellCenterX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
        RDouble4D &gradUVWTCellCenterY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
        RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
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
                    RDouble S = two * (s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23));
                    S = sqrt(S);
                    (*strainRate)(i, j, k) = S;

                    RDouble Omx = dwdy - dvdz;
                    RDouble Omy = dudz - dwdx;
                    RDouble Omz = dvdx - dudy;
                    RDouble Om = DISTANCE(Omx, Omy, Omz);

                    (*XVorticity)(i, j, k) = Omx;
                    (*YVorticity)(i, j, k) = Omy;
                    (*ZVorticity)(i, j, k) = Omz;
                    (*vorticityMagnitude)(i, j, k) = Om;

                    (*QCriteria)(i, j, k) = fourth * (Om * Om - S * S);

                }
            }
        }
    }

    void NSSolverStruct::ComputeCp(Grid *gridIn, RDouble4D *&pressureCoef)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Range I, J, K;
        GetRange(ni, nj, nk, -2, 1, I, J, K);

        Range M(0, 0);

        pressureCoef = new RDouble4D(I, J, K, M, fortranArray);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        RDouble4D &primitiveVars = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        using namespace IDX;
        Param_NSSolverStruct *parameters = GetControlParameters();
        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
        RDouble presureFarfield = primitiveVarFarfield[IP];

        int m = 0;
        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    (*pressureCoef)(i, j, k, m) = two * (primitiveVars(i, j, k, IP) - presureFarfield);
                }
            }
        }
    }

    RDouble4D *NSSolverStruct::ComputeDimensionalVariables(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Range I, J, K;
        GetRange(ni, nj, nk, -2, 1, I, J, K);

        //Allocate the memories for density, pressure, temperature and the number density of electron.
        Range M(0, 5);
        RDouble4D *dimensionalVariables = new RDouble4D(I, J, K, M, fortranArray);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        RDouble4D &primitiveVars = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));

        using namespace IDX;
        using namespace GAS_SPACE;
        Param_NSSolverStruct *parameters = GetControlParameters();
        RDouble refDensity = parameters->GetRefDimensionalDensity();
        RDouble refTemperature = parameters->GetRefDimensionalTemperature();
        RDouble refVelocity = parameters->GetRefDimensionalVelocity();//parameters->GetRefMachNumber() * parameters->GetRefDimensionalSonicSpeed();
        RDouble refPressure = refDensity * refVelocity * refVelocity;

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    //! Compute the dimensional density.
                    (*dimensionalVariables)(i, j, k, 0) = refDensity * primitiveVars(i, j, k, IR);

                    //! Compute the dimensional velocity.
                    (*dimensionalVariables)(i, j, k, 1) = refVelocity * primitiveVars(i, j, k, IU);
                    (*dimensionalVariables)(i, j, k, 2) = refVelocity * primitiveVars(i, j, k, IV);
                    (*dimensionalVariables)(i, j, k, 3) = refVelocity * primitiveVars(i, j, k, IW);

                    //! Compute the dimensional pressure.
                    (*dimensionalVariables)(i, j, k, 4) = refPressure * primitiveVars(i, j, k, IP);

                    //! Compute the dimensional temperature.
                    (*dimensionalVariables)(i, j, k, 5) = refTemperature * temperatures(i, j, k, ITT);
                }
            }
        }

        return dimensionalVariables;
    }

    RDouble4D *NSSolverStruct::ComputeDimensionalElectron(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Range I, J, K;
        GetRange(ni, nj, nk, -2, 1, I, J, K);

        //Allocate the memories for density, pressure, temperature and the number density of electron.
        Range M(0, 0);
        RDouble4D *tmpArray = new RDouble4D(I, J, K, M, fortranArray);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        RDouble4D &primitiveVars = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        using namespace IDX;
        using namespace GAS_SPACE;
        Param_NSSolverStruct *parameters = GetControlParameters();
        RDouble refDensity = parameters->GetRefDimensionalDensity();
        int nNSEquation = parameters->GetNSEquationNumber();
        int indexOfElectron = parameters->GetIndexOfElectron();

        RDouble rho = 0.0, Ce = 0.0;
        if (indexOfElectron >= 0)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        //Compute the dimensional density.
                        rho = refDensity * primitiveVars(i, j, k, IR);

                        //Compute the number density of electron.
                        Ce = primitiveVars(i, j, k, nNSEquation + indexOfElectron);

                        (*tmpArray)(i, j, k, 0) = MAX(1.0e-25, rho * Ce * AVOGADRO_CONSTANT / 5.486e-7);
                    }
                }
            }
        }
        else
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        (*tmpArray)(i, j, k, 0) = 0.0;
                    }
                }
            }
        }

        return tmpArray;
    }

    RDouble4D *NSSolverStruct::ComputeMachNumberField(Grid *gridIn)
    {
        const RDouble gamma = 1.4;
        StructGrid *grid = StructGridCast(gridIn);

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Range I, J, K;
        GetRange(ni, nj, nk, -2, 1, I, J, K);
        Range M(0, 0);

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        RDouble4D *mach = new RDouble4D(I, J, K, M, fortranArray);

        using namespace IDX;

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        int m = 0;
        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    RDouble rm = q(i, j, k, IR);
                    RDouble um = q(i, j, k, IU);
                    RDouble vm = q(i, j, k, IV);
                    RDouble wm = q(i, j, k, IW);
                    RDouble pm = q(i, j, k, IP);
                    RDouble c2 = abs(gamma * pm / (rm + SMALL));
                    RDouble v2 = um * um + vm * vm + wm * wm;

                    (*mach)(i, j, k, m) = sqrt(v2 / (c2 + SMALL));
                }
            }
        }

        return mach;
    }

    RDouble4D *NSSolverStruct::CompMachNumber(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Range I, J, K;
        GetRange(ni, nj, nk, -2, 1, I, J, K);
        Range M(0, 0);

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble3D &gama = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));
        RDouble4D *mach = new RDouble4D(I, J, K, M, fortranArray);

        using namespace IDX;

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        int m = 0;
        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    RDouble rm = q(i, j, k, IR);
                    RDouble um = q(i, j, k, IU);
                    RDouble vm = q(i, j, k, IV);
                    RDouble wm = q(i, j, k, IW);
                    RDouble pm = q(i, j, k, IP);
                    RDouble c2 = abs(gama(i, j, k) * pm / (rm + SMALL));
                    RDouble v2 = um * um + vm * vm + wm * wm;

                    (*mach)(i, j, k, m) = sqrt(v2 / (c2 + SMALL));
                }
            }
        }

        return mach;
    }

    void NSSolverStruct::VisualizationAverageFlow(ActionKey *actkey)
    {

    }

    //! Restrict from fine to coarse grid for all q.
    void NSSolverStruct::RestrictAllQ(Grid *fineGridIn, Grid *coarseGridIn)
    {
        StructGrid *fineGrid = StructGridCast(fineGridIn);
        StructGrid *coarseGrid = StructGridCast(coarseGridIn);

        int coarseNi = coarseGrid->GetNI();
        int coarseNj = coarseGrid->GetNJ();
        int coarseNk = coarseGrid->GetNK();

        RDouble4D &finePrimitiveVars = *reinterpret_cast <RDouble4D *> (fineGrid->GetDataPtr("q"));
        RDouble4D &coarsePrimitiveVars = *reinterpret_cast <RDouble4D *> (coarseGrid->GetDataPtr("q"));

        Param_NSSolverStruct *parameters = GetControlParameters();

        int nEquation = GetNumberOfEquations();

        //! Zero Setting.

        coarsePrimitiveVars = 0.0;
        RestrictQ(fineGrid, finePrimitiveVars, coarseGrid, coarsePrimitiveVars, nEquation);
        //coarseGrid->GhostCell3DExceptInterface(coarsePrimitiveVars, nEquation);
        //GhostCell3D(coarsePrimitiveVars, coarseNi, coarseNj, coarseNk, nEquation);

        RDouble3D &fineGama = *reinterpret_cast <RDouble3D *> (fineGrid->GetDataPtr("gama"));
        RDouble3D &coarseGama = *reinterpret_cast <RDouble3D *> (coarseGrid->GetDataPtr("gama"));
        RestrictQ(fineGrid, fineGama, coarseGrid, coarseGama);
        //GhostCell3D(coarseGama, coarseNi, coarseNj, coarseNk);

        int viscousType = parameters->GetViscousType();
        if (viscousType > INVISCID)
        {
            RDouble3D &fineViscosityLaminar = *reinterpret_cast <RDouble3D *> (fineGrid->GetDataPtr("visl"));
            RDouble3D &coarseViscosityLaminar = *reinterpret_cast <RDouble3D *> (coarseGrid->GetDataPtr("visl"));
            coarseViscosityLaminar = 0.0;
            RestrictQ(fineGrid, fineViscosityLaminar, coarseGrid, coarseViscosityLaminar);

            GhostCell3D(coarseViscosityLaminar, coarseNi, coarseNj, coarseNk);

            if (viscousType > LAMINAR)
            {
                RDouble3D &fineViscosityTurbulence = *reinterpret_cast <RDouble3D *> (fineGrid->GetDataPtr("vist"));
                RDouble3D &coarseViscosityTurbulence = *reinterpret_cast <RDouble3D *> (coarseGrid->GetDataPtr("vist"));
                coarseViscosityTurbulence = 0.0;
                RestrictQ(fineGrid, fineViscosityTurbulence, coarseGrid, coarseViscosityTurbulence);

                GhostCell3D(coarseViscosityTurbulence, coarseNi, coarseNj, coarseNk);
                //coarseGrid->GhostCell3DExceptInterface(coarseViscosityTurbulence);
            }
        }
    }

    //! Restrict defect from fine to coarse grid.
    void NSSolverStruct::RestrictDefect(Grid *fineGridIn, Grid *coarseGridIn)
    {
        StructGrid *fineGrid = StructGridCast(fineGridIn);
        StructGrid *coarseGrid = StructGridCast(coarseGridIn);

        //! The following statement is important.
        ZeroResiduals(coarseGrid);

        RDouble4D &coarseResidual = *reinterpret_cast <RDouble4D *> (coarseGrid->GetDataPtr("res"));
        RDouble4D &fineResidual = *reinterpret_cast <RDouble4D *> (fineGrid->GetDataPtr("res"));

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nLaminar = parameters->GetLaminarNumber();
        int nEquation = GetNumberOfEquations();

        int coarseNi = coarseGrid->GetNI();
        int coarseNj = coarseGrid->GetNJ();
        int coarseNk = coarseGrid->GetNK();

        int iMultigrid = fineGrid->GetMultigridStepI();
        int jMultigrid = fineGrid->GetMultigridStepJ();
        int kMultigrid = fineGrid->GetMultigridStepK();

        int id = iMultigrid / 2;
        int jd = jMultigrid / 2;
        int kd = kMultigrid / 2;

        Range I(1, MAX(coarseNi - 1, 1));
        Range J(1, MAX(coarseNj - 1, 1));
        Range K(1, MAX(coarseNk - 1, 1));
        Range M(0, nEquation - 1);

        RDouble coeff = eighth * static_cast<RDouble>(iMultigrid * jMultigrid * kMultigrid);

        int i0, j0, k0, i, j, k, ip, jp, kp;
        for (k0 = K.first(); k0 <= K.last(); ++ k0)
        {
            k = (k0 - 1) * kMultigrid + 1;
            kp = k + kd;
            for (j0 = J.first(); j0 <= J.last(); ++ j0)
            {
                j = (j0 - 1) * jMultigrid + 1;
                jp = j + jd;
                for (i0 = I.first(); i0 <= I.last(); ++ i0)
                {
                    i = (i0 - 1) * iMultigrid + 1;
                    ip = i + id;
                    for (int iLaminar = 0; iLaminar < nLaminar; ++ iLaminar)
                    {
                        coarseResidual(i0, j0, k0, iLaminar) -= coeff * (fineResidual(i, j, k, iLaminar) + fineResidual(ip, j, k, iLaminar)
                            + fineResidual(ip, jp, k, iLaminar) + fineResidual(i, jp, k, iLaminar)
                            + fineResidual(i, j, kp, iLaminar) + fineResidual(ip, j, kp, iLaminar)
                            + fineResidual(ip, jp, kp, iLaminar) + fineResidual(i, jp, kp, iLaminar));
                    }
                }
            }
        }
    }

    //! Put correction back on the coarse grid.
    void NSSolverStruct::PutCorrectionBack(Grid *grid, FieldProxy *qProxy)
    {
        RDouble4D &qold = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &q = qProxy->GetField_STR();
        qold += q;
    }

    //! Correct variable in fine grid using correction in coarse grid.
    void NSSolverStruct::CorrectFineGrid(Grid *fineGridIn, Grid *coarseGridIn)
    {
        StructGrid *fineGrid = StructGridCast(fineGridIn);
        StructGrid *coarseGrid = StructGridCast(coarseGridIn);

        RDouble4D &finePrimitiveVars = *reinterpret_cast <RDouble4D *> (fineGrid->GetDataPtr("q"));

        int fineNi = fineGrid->GetNI();
        int fineNj = fineGrid->GetNJ();
        int fineNk = fineGrid->GetNK();

        Param_NSSolverStruct *parameters = GetControlParameters();
        int mgProlongationType = parameters->GetMgProlongationType();

        int nEquation = GetNumberOfEquations();

        Range I, J, K;
        GetRange(fineNi, fineNj, fineNk, -2, 1, I, J, K);
        Range M(0, nEquation - 1);

        if (fineNk == 1)
        {
            K.setRange(1, 1);
        }

        RDouble4D work(I, J, K, M, fortranArray);

        FieldProxy *workProxy = new FieldProxy();
        workProxy->SetField_STR(&work);

        FieldProxy *finePrimitiveVarsProxy = new FieldProxy();
        finePrimitiveVarsProxy->SetField_STR(&finePrimitiveVars);

        if (mgProlongationType == 1)
        {
            prol(fineGrid, coarseGrid, workProxy);
        }
        else if (mgProlongationType == 2)
        {
            prolhic_eric(fineGrid, coarseGrid, workProxy);
        }

        //! Perform 'w = w + work' to update the state vector.
        PositiveUpdate(fineGrid, finePrimitiveVarsProxy, workProxy);

        delete finePrimitiveVarsProxy;
        delete workProxy;
    }

    void NSSolverStruct::InterpolatFineGrid(Grid *fineGridIn, Grid *coarseGridIn)
    {
        StructGrid *fineGrid = StructGridCast(fineGridIn);
        StructGrid *coarseGrid = StructGridCast(coarseGridIn);

        RDouble4D &finePrimitiveVars = *reinterpret_cast <RDouble4D *> (fineGrid->GetDataPtr("q"));
        RDouble4D &coarsePrimitiveVars = *reinterpret_cast <RDouble4D *> (coarseGrid->GetDataPtr("q"));

        Param_NSSolverStruct *parameters = GetControlParameters();

        int nEquation = GetNumberOfEquations();

        InterpolatQ(fineGrid, finePrimitiveVars, coarseGrid, coarsePrimitiveVars, nEquation);

        RDouble3D &fineGama = *reinterpret_cast <RDouble3D *> (fineGrid->GetDataPtr("gama"));
        RDouble3D &coarseGama = *reinterpret_cast <RDouble3D *> (coarseGrid->GetDataPtr("gama"));
        InterpolatQ(fineGrid, fineGama, coarseGrid, coarseGama);

        int viscousType = parameters->GetViscousType();
        if (viscousType > INVISCID)
        {
            RDouble3D &fineViscosityLaminar = *reinterpret_cast <RDouble3D *> (fineGrid->GetDataPtr("visl"));
            RDouble3D &coarseViscosityLaminar = *reinterpret_cast <RDouble3D *> (coarseGrid->GetDataPtr("visl"));
            InterpolatQ(fineGrid, fineViscosityLaminar, coarseGrid, coarseViscosityLaminar);

            if (viscousType > LAMINAR)
            {
                RDouble3D &fineViscosityTurbulence = *reinterpret_cast <RDouble3D *> (fineGrid->GetDataPtr("vist"));
                RDouble3D &coarseViscosityTurbulence = *reinterpret_cast <RDouble3D *> (coarseGrid->GetDataPtr("vist"));
                InterpolatQ(fineGrid, fineViscosityTurbulence, coarseGrid, coarseViscosityTurbulence);
            }
        }
    }

    void NSSolverStruct::prol(Grid *fineGridIn, Grid *coarseGridIn, FieldProxy *workProxy)
    {
        StructGrid *fineGrid = StructGridCast(fineGridIn);
        StructGrid *coarseGrid = StructGridCast(coarseGridIn);

        RDouble4D &work = workProxy->GetField_STR();

        RDouble4D &coarsePrimitiveVars = *reinterpret_cast <RDouble4D *> (coarseGrid->GetDataPtr("q"));

        int nEquation = GetNumberOfEquations();

        int fineNi = fineGrid->GetNI();
        int fineNj = fineGrid->GetNJ();
        int fineNk = fineGrid->GetNK();

        int coarseNi = coarseGrid->GetNI();
        int coarseNj = coarseGrid->GetNJ();
        int coarseNk = coarseGrid->GetNK();

        int iMultigridStep = fineGrid->GetMultigridStepI();
        int jMultigridStep = fineGrid->GetMultigridStepJ();
        int kMultigridStep = fineGrid->GetMultigridStepK();

        int id = iMultigridStep / 2;
        int jd = jMultigridStep / 2;
        int kd = kMultigridStep / 2;

        Range I(1, MAX(coarseNi - 1, 1));
        Range J(1, MAX(coarseNj - 1, 1));
        Range K(1, MAX(coarseNk - 1, 1));

        int i0, j0, k0, i, j, k, ip, jp, kp;
        for (k0 = K.first(); k0 <= K.last(); ++ k0)
        {
            k = (k0 - 1) * kMultigridStep + 1;
            kp = k + kd;
            for (j0 = J.first(); j0 <= J.last(); ++ j0)
            {
                j = (j0 - 1) * jMultigridStep + 1;
                jp = j + jd;
                for (i0 = I.first(); i0 <= I.last(); ++ i0)
                {
                    i = (i0 - 1) * iMultigridStep + 1;
                    ip = i + id;
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        RDouble w0diff = coarsePrimitiveVars(i0, j0, k0, m);

                        work(i, j, k, m) = w0diff;
                        work(ip, j, k, m) = w0diff;
                        work(ip, jp, k, m) = w0diff;
                        work(i, jp, k, m) = w0diff;
                        work(i, j, kp, m) = w0diff;
                        work(ip, j, kp, m) = w0diff;
                        work(ip, jp, kp, m) = w0diff;
                        work(i, jp, kp, m) = w0diff;
                    }
                }
            }
        }

        GhostCell3D(work, fineNi, fineNj, fineNk, nEquation);
    }

    void NSSolverStruct::prolhic_eric(Grid *fineGridIn, Grid *coarseGridIn, FieldProxy *workProxy)
    {
        StructGrid *fineGrid = StructGridCast(fineGridIn);
        StructGrid *coarseGrid = StructGridCast(coarseGridIn);

        RDouble4D &work = workProxy->GetField_STR();
        RDouble4D &coarsePrimitiveVars = *reinterpret_cast <RDouble4D *> (coarseGrid->GetDataPtr("q"));

        int nEquation = GetNumberOfEquations();

        int fineNi = fineGrid->GetNI();
        int fineNj = fineGrid->GetNJ();
        int fineNk = fineGrid->GetNK();

        int coarseNi = coarseGrid->GetNI();
        int coarseNj = coarseGrid->GetNJ();
        int coarseNk = coarseGrid->GetNK();

        int iMultigridStep = fineGrid->GetMultigridStepI();
        int jMultigridStep = fineGrid->GetMultigridStepJ();
        int kMultigridStep = fineGrid->GetMultigridStepK();

        int id = iMultigridStep / 2;
        int jd = jMultigridStep / 2;
        int kd = kMultigridStep / 2;

        int nDim = fineGrid->GetDim();

        Range I, J, K;
        GetRange(fineNi, fineNj, fineNk, -2, 1, I, J, K);

        Range CoarseI, CoarseJ, CoarseK;
        GetRange(coarseNi, coarseNj, coarseNk, -2, 1, CoarseI, CoarseJ, CoarseK);

        RDouble3D wq(CoarseI, CoarseJ, CoarseK, fortranArray);

        RDouble3D wdiff1(I, J, K, fortranArray);
        RDouble3D wdiff2(I, J, K, fortranArray);

        const RDouble v0 = 0.25;
        const RDouble v1 = 0.75;

        if (nDim == TWO_D)
        {
            //! Two-dimensional computations.

            int ist, ied, jst, jed, kst, ked;
            ist = 1;
            ied = coarseNi - 1;
            jst = 1;
            jed = coarseNj - 1;
            kst = 1;
            ked = 1;

            for (int m = 0; m < nEquation; ++ m)
            {
                for (int kk = kst; kk <= ked; ++ kk)
                {
                    for (int jj = jst; jj <= jed; ++ jj)
                    {
                        for (int ii = ist; ii <= ied; ++ ii)
                        {
                            wq(ii, jj, kk) = coarsePrimitiveVars(ii, jj, kk, m);
                        }
                    }
                }

                GhostCell3D(wq, coarseNi, coarseNj, coarseNk);
                FillCornerPoint3D(wq, coarseNi, coarseNj, coarseNk);

                //! Interpolation in I direction.
                int k0 = 1;
                int j = 0;

                for (int j0 = 0; j0 <= coarseNj; ++ j0)
                {
                    int i = 0;
                    for (int i0 = 0; i0 <= coarseNi - 1; ++ i0)
                    {
                        RDouble w0diff0 = wq(i0, j0, k0);
                        RDouble w0diff1 = wq(i0 + id, j0, k0);

                        wdiff1(i, j0, k0) = v1 * w0diff0 + v0 * w0diff1;
                        wdiff1(i + id, j0, k0) = v0 * w0diff0 + v1 * w0diff1;

                        i += 1 + id;
                    }
                    j += 1 + jd;
                }

                //! Interpolation in J direction.
                for (int i = 0; i <= fineNi; ++ i)
                {
                    j = 0;
                    for (int j0 = 0; j0 <= coarseNj - 1; ++ j0)
                    {
                        RDouble w0diff0 = wdiff1(i, j0, k0);
                        RDouble w0diff1 = wdiff1(i, j0 + jd, k0);
                        work(i, j, k0, m) = v1 * w0diff0 + v0 * w0diff1;
                        work(i, j + jd, k0, m) = v0 * w0diff0 + v1 * w0diff1;
                        j += 1 + jd;
                    }
                }

            }
        }
        else
        {

            int ist, ied, jst, jed, kst, ked;
            ist = 1;
            ied = coarseNi - 1;
            jst = 1;
            jed = coarseNj - 1;
            kst = 1;
            ked = coarseNk - 1;

            work = 0.0;

            for (int m = 0; m < nEquation; ++ m)
            {
                for (int kk = kst; kk <= ked; ++ kk)
                {
                    for (int jj = jst; jj <= jed; ++ jj)
                    {
                        for (int ii = ist; ii <= ied; ++ ii)
                        {
                            wq(ii, jj, kk) = coarsePrimitiveVars(ii, jj, kk, m);
                        }
                    }
                }

                GhostCell3D(wq, coarseNi, coarseNj, coarseNk);
                FillCornerPoint3D(wq, coarseNi, coarseNj, coarseNk);

                //! Interpolation in I direction.
                for (int k0 = 0; k0 <= coarseNk; ++ k0)
                {
                    for (int j0 = 0; j0 <= coarseNj; ++ j0)
                    {
                        int i = 0;
                        for (int i0 = 0; i0 <= coarseNi - 1; ++ i0)
                        {
                            RDouble w0diff0 = wq(i0, j0, k0);
                            RDouble w0diff1 = wq(i0 + id, j0, k0);

                            wdiff1(i, j0, k0) = v1 * w0diff0 + v0 * w0diff1;
                            wdiff1(i + id, j0, k0) = v0 * w0diff0 + v1 * w0diff1;

                            work(i, j0, k0, m) = v1 * w0diff0 + v0 * w0diff1;
                            work(i + id, j0, k0, m) = v0 * w0diff0 + v1 * w0diff1;

                            i += 1 + id;
                        }
                    }
                }

                //! Interpolation in J direction.
                for (int k0 = 0; k0 <= coarseNk; ++ k0)
                {
                    for (int i = 0; i <= fineNi; ++ i)
                    {
                        int j = 0;
                        for (int j0 = 0; j0 <= coarseNj - 1; ++ j0)
                        {
                            RDouble w0diff0 = wdiff1(i, j0, k0);
                            RDouble w0diff1 = wdiff1(i, j0 + jd, k0);

                            wdiff2(i, j, k0) = v1 * w0diff0 + v0 * w0diff1;
                            wdiff2(i, j + jd, k0) = v0 * w0diff0 + v1 * w0diff1;

                            work(i, j, k0, m) = v1 * w0diff0 + v0 * w0diff1;
                            work(i, j + jd, k0, m) = v0 * w0diff0 + v1 * w0diff1;
                            j += 1 + jd;
                        }
                    }
                }

                //! Interpolation in K direction.
                for (int j = 0; j <= fineNj; ++ j)
                {
                    for (int i = 0; i <= fineNi; ++ i)
                    {
                        int k = 0;
                        for (int k0 = 0; k0 <= coarseNk - 1; ++ k0)
                        {
                            RDouble w0diff0 = wdiff2(i, j, k0);
                            RDouble w0diff1 = wdiff2(i, j, k0 + kd);

                            work(i, j, k, m) = v1 * w0diff0 + v0 * w0diff1;
                            work(i, j, k + kd, m) = v0 * w0diff0 + v1 * w0diff1;

                            k += 1 + kd;
                        }
                    }
                }
            }
        }
    }

    void NSSolverStruct::PositiveUpdate(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &q = qProxy->GetField_STR();
        RDouble4D &dq = dqProxy->GetField_STR();

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        RDouble alpha = half;

        using namespace IDX;

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    RDouble rm = q(i, j, k, IR);
                    RDouble um = q(i, j, k, IU);
                    RDouble vm = q(i, j, k, IV);
                    RDouble wm = q(i, j, k, IW);
                    RDouble pm = q(i, j, k, IP);

                    rm += alpha * dq(i, j, k, IR);
                    um += alpha * dq(i, j, k, IU);
                    vm += alpha * dq(i, j, k, IV);
                    wm += alpha * dq(i, j, k, IW);
                    pm += alpha * dq(i, j, k, IP);

                    if (pm <= 0.0 || rm <= 0.0)
                    {
                        //! Let dt be one order smaller.
                        RDouble coef = 0.9;
                        rm -= dq(i, j, k, IR) * coef;
                        um -= dq(i, j, k, IU) * coef;
                        vm -= dq(i, j, k, IV) * coef;
                        wm -= dq(i, j, k, IW) * coef;
                        pm -= dq(i, j, k, IP) * coef;

                        if (pm <= 0.0 || rm <= 0.0)
                        {
                            coef = 0.09;
                            rm -= dq(i, j, k, IR) * coef;
                            um -= dq(i, j, k, IU) * coef;
                            vm -= dq(i, j, k, IV) * coef;
                            wm -= dq(i, j, k, IW) * coef;
                            pm -= dq(i, j, k, IP) * coef;
                        }
                    }

                    if (pm > 0.0 && rm > 0.0)
                    {
                        q(i, j, k, IR) = rm;
                        q(i, j, k, IU) = um;
                        q(i, j, k, IV) = vm;
                        q(i, j, k, IW) = wm;
                        q(i, j, k, IP) = pm;
                    }
                }
            }
        }
    }

    void NSSolverStruct::VisualizationVariables(Grid *gridIn, FieldProxy *workProxy, const string &filename)
    {
        StructGrid *grid = StructGridCast(gridIn);

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int nEquation = GetNumberOfEquations();

        vector <string> title_tecplot;
        title_tecplot.push_back("title=\"Flow Fields of PHengLEI\"");
        title_tecplot.push_back("variables=");
        title_tecplot.push_back("\"x\"");
        title_tecplot.push_back("\"y\"");
        title_tecplot.push_back("\"z\"");
        title_tecplot.push_back("\"work1\"");
        title_tecplot.push_back("\"work2\"");
        title_tecplot.push_back("\"work3\"");
        title_tecplot.push_back("\"work4\"");
        title_tecplot.push_back("\"work5\"");

        int nvarplot;
        nvarplot = 10;

        Range II, JJ, KK;
        GetRange(ni, nj, nk, 0, 0, II, JJ, KK);
        Range MM(0, nEquation - 1);

        RDouble4D workn(II, JJ, KK, MM, fortranArray);

        FieldProxy *worknProxy = new FieldProxy();

        worknProxy->SetField_STR(&workn);
        CompNodeVar(grid, worknProxy, workProxy);

        delete worknProxy;

        RDouble3D &xx = *grid->GetStructX();
        RDouble3D &yy = *grid->GetStructY();
        RDouble3D &zz = *grid->GetStructZ();

        using namespace PHSPACE;

        std::ostringstream oss;

        for (std::size_t i = 0; i < title_tecplot.size(); ++ i)
        {
            oss << title_tecplot[i] << "\n";
        }

        int nTotalCell = grid->GetNTotalCell();
        int nTotalNode = grid->GetNTotalNode();

        int wordwidth = 20;

        if (nk == 1)
        {
            oss << "zone  N = " << nTotalNode << "  E = " << nTotalCell << " f = FEPOINT, ET = QUADRILATERAL\n";

            int iNodeStart, iNodeEnd, jNodeStart, jNodeEnd, kNodeStart, kNodeEnd;
            grid->GetNodeIterationIndex(iNodeStart, iNodeEnd, jNodeStart, jNodeEnd, kNodeStart, kNodeEnd);

            for (int k = kNodeStart; k <= kNodeEnd; ++ k)
            {
                for (int j = jNodeStart; j <= jNodeEnd; ++ j)
                {
                    for (int i = iNodeStart; i <= iNodeEnd; ++ i)
                    {
                        oss << setiosflags(ios::left);
                        oss << setiosflags(ios::scientific);
                        oss << setprecision(10);
                        oss << setw(wordwidth) << xx(i, j, k)
                            << setw(wordwidth) << yy(i, j, k)
                            << setw(wordwidth) << zz(i, j, k);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            oss << setw(wordwidth) << workn(i, j, k, m);
                        }
                        oss << "\n";
                    }
                }
            }

            for (int j = 1; j <= nj - 1; ++ j)
            {
                for (int i = 1; i <= ni - 1; ++ i)
                {
                    oss << i + (j - 1) * ni << " ";
                    oss << i + 1 + (j - 1) * ni << " ";
                    oss << i + 1 + (j)*ni << " ";
                    oss << i + (j)*ni << "\n";
                }
            }
        }

        int s_st[3], s_ed[3];

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();
            bool flag;
            flag = (IsWall(BCType))
                || (BCType == PHENGLEI::SYMMETRY)
                || (BCType == PHENGLEI::OUTFLOW);

            if (!flag)
            {
                continue;
            }

            int *s_lr3d = structBC->GetFaceDirectionIndex();
            int nDim = GetDim();

            for (int m = 0; m < nDim; ++ m)
            {
                s_st[m] = structBC->GetStartPoint(m);
                s_ed[m] = structBC->GetEndPoint(m);
                if (structBC->GetFaceDirection() == m)
                {
                    if (s_st[m] > 1 || s_lr3d[m] == 1)
                    {
                        s_st[m] += 1;
                        s_ed[m] += 1;
                    }
                }
                else
                {
                    s_ed[m] += 1;
                }
            }

            int ist = s_st[0];
            int ied = s_ed[0];
            int jst = s_st[1];
            int jed = s_ed[1];
            int kst = s_st[2];
            int ked = s_ed[2];

            Range I, J, K;
            GetRange(I, J, K, ist, ied, jst, jed, kst, ked);

            if (nDim != THREE_D)
            {
                K.setRange(1, 1);
                kst = 1;
                ked = 1;
            }

            Int3D newindex(I, J, K, fortranArray);

            int nSurface = structBC->GetFaceDirection() + 1;

            int count = 0;
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        newindex(i, j, k) = count;
                        ++ count;
                    }
                }
            }

            int nTotalNode = count;
            int nTotalCell = 1;

            for (int m = 0; m < 3; ++ m)
            {
                nTotalCell *= (structBC->GetEndPoint(m) - structBC->GetStartPoint(m) + 1);
            }

            oss << "zone N = " << nTotalNode << "  E = " << nTotalCell << " f = FEPOINT, ET = QUADRILATERAL\n";

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        oss << setiosflags(ios::left);
                        oss << setiosflags(ios::scientific);
                        oss << setprecision(10);
                        oss << setw(wordwidth) << xx(i, j, k)
                            << setw(wordwidth) << yy(i, j, k)
                            << setw(wordwidth) << zz(i, j, k);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            oss << setw(wordwidth) << workn(i, j, k, m);
                        }
                        oss << "\n";
                    }
                }
            }

            int il1, jl1, kl1;
            grid->GetNsurfIndex(il1, jl1, kl1, nSurface);


            il1 = 1 - il1;
            jl1 = 1 - jl1;
            kl1 = 1 - kl1;

            if (nk == 1)
            {
                kl1 = 0;
            }

            for (int k = kst; k <= ked - kl1; ++ k)
            {
                for (int j = jst; j <= jed - jl1; ++ j)
                {
                    for (int i = ist; i <= ied - il1; ++ i)
                    {
                        int p1, p2, p3, p4;

                        int il, jl, kl;
                        grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, nSurface);

                        if (nSurface == 1)
                        {
                            p1 = newindex(i, j, k) + 1;
                            p2 = newindex(i, jl, k) + 1;
                            p3 = newindex(i, jl, kl) + 1;
                            p4 = newindex(i, j, kl) + 1;
                        }
                        else if (nSurface == 2)
                        {
                            p1 = newindex(i, j, k) + 1;
                            p2 = newindex(il, j, k) + 1;
                            p3 = newindex(il, j, kl) + 1;
                            p4 = newindex(i, j, kl) + 1;
                        }
                        else
                        {
                            p1 = newindex(i, j, k) + 1;
                            p2 = newindex(il, j, k) + 1;
                            p3 = newindex(il, jl, k) + 1;
                            p4 = newindex(i, jl, k) + 1;
                        }
                        oss << p1 << " ";
                        oss << p2 << " ";
                        oss << p3 << " ";
                        oss << p4 << "\n";
                    }
                }
            }
        }

        fstream file;
        file.open(filename.c_str(), ios_base::out | ios_base::trunc);
        if (!file)
        {
            TK_Exit::FileOpenErrorExit(filename);
        }

        file << oss.str();
        file.close();
    }

    void NSSolverStruct::CompNodeVar(Grid *gridIn, RDouble4D &qn, int m, RDouble4D &q, int n)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int nk = grid->GetNK();

        int il1 = 1;
        int jl1 = 1;
        int kl1 = 1;
        if (nk == 1)
        {
            kl1 = 0;
        }

        int il, jl, kl;

        int iNodeStart, iNodeEnd, jNodeStart, jNodeEnd, kNodeStart, kNodeEnd;
        grid->GetNodeIterationIndex(iNodeStart, iNodeEnd, jNodeStart, jNodeEnd, kNodeStart, kNodeEnd);

        for (int k = kNodeStart; k <= kNodeEnd; ++ k)
        {
            for (int j = jNodeStart; j <= jNodeEnd; ++ j)
            {
                for (int i = iNodeStart; i <= iNodeEnd; ++ i)
                {
                    il = i - il1;
                    jl = j - jl1;
                    kl = k - kl1;
                    qn(i, j, k, m) = eighth * (q(il, jl, kl, n) + q(il, j, kl, n) + q(i, jl, kl, n) + q(i, j, kl, n) +
                        q(il, jl, k, n) + q(il, j, k, n) + q(i, jl, k, n) + q(i, j, k, n));
                }
            }
        }
    }

    void NSSolverStruct::CompNodeVar(Grid *gridIn, RDouble4D &qn, int m, RDouble3D &q)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int nk = grid->GetNK();

        int il1 = 1;
        int jl1 = 1;
        int kl1 = 1;
        if (nk == 1)
        {
            kl1 = 0;
        }

        int il, jl, kl;

        int iNodeStart, iNodeEnd, jNodeStart, jNodeEnd, kNodeStart, kNodeEnd;
        grid->GetNodeIterationIndex(iNodeStart, iNodeEnd, jNodeStart, jNodeEnd, kNodeStart, kNodeEnd);

        for (int k = kNodeStart; k <= kNodeEnd; ++ k)
        {
            for (int j = jNodeStart; j <= jNodeEnd; ++ j)
            {
                for (int i = iNodeStart; i <= iNodeEnd; ++ i)
                {
                    il = i - il1;
                    jl = j - jl1;
                    kl = k - kl1;
                    qn(i, j, k, m) = eighth * (q(il, jl, kl) + q(il, j, kl) + q(i, jl, kl) + q(i, j, kl) +
                        q(il, jl, k) + q(il, j, k) + q(i, jl, k) + q(i, j, k));
                }
            }
        }
    }

    void NSSolverStruct::CompNodeVar(Grid *gridIn, FieldProxy *qnProxy, FieldProxy *qProxy)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int nk = grid->GetNK();

        RDouble4D &q = qProxy->GetField_STR();
        RDouble4D &q_n = qnProxy->GetField_STR();

        int nEquation = GetNumberOfEquations();

        int il1 = 1;
        int jl1 = 1;
        int kl1 = 1;
        if (nk == 1)
        {
            kl1 = 0;
        }

        int il, jl, kl;

        int iNodeStart, iNodeEnd, jNodeStart, jNodeEnd, kNodeStart, kNodeEnd;
        grid->GetNodeIterationIndex(iNodeStart, iNodeEnd, jNodeStart, jNodeEnd, kNodeStart, kNodeEnd);

        Range M(0, nEquation - 1);
        int mst = M.first();
        int med = M.last();

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kNodeStart; k <= kNodeEnd; ++ k)
            {
                for (int j = jNodeStart; j <= jNodeEnd; ++ j)
                {
                    for (int i = iNodeStart; i <= iNodeEnd; ++ i)
                    {
                        il = i - il1;
                        jl = j - jl1;
                        kl = k - kl1;
                        q_n(i, j, k, m) = eighth * (q(il, jl, kl, m) + q(il, j, kl, m) + q(i, jl, kl, m) + q(i, j, kl, m) +
                            q(il, jl, k, m) + q(il, j, k, m) + q(i, jl, k, m) + q(i, j, k, m));
                    }
                }
            }
        }

        InterfaceValueFix(gridIn, qnProxy, qProxy);
        WallValueFix(gridIn, qnProxy, qProxy);
    }

    void NSSolverStruct::InterfaceValueFix(Grid *gridIn, FieldProxy *qnProxy, FieldProxy *qProxy)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &q = qProxy->GetField_STR();
        RDouble4D &q_n = qnProxy->GetField_STR();

        int nEquation = GetNumberOfEquations();

        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        if (!interfaceInfo)
        {
            return;
        }

        int nTotalNode = grid->GetNTotalNode();

        int nlist = 4;
        if (grid->GetDim() == TWO_D)
        {
            nlist = 2;
        }

        int ijk[4][3];
        int irot[3], jrot[3], krot[3];

        GetIJKROT(ijk, irot, jrot, krot);

        //VVInt node_list;
        //set< PHSort< VInt > > node_sort;
        //VInt  nodeindex(3);

        int ist, ied, jst, jed, kst, ked;
        int iSurface, BCType, index;
        int il, jl, kl;
        int is, js, ks;
        int it, jt, kt;
        int ir, jr, kr;
        int in, jn, kn;

        RDouble *qtmp = new RDouble[nEquation];
        int *n_count = new int[nTotalNode];
        RDouble **q_n_list = NewPointer2<RDouble>(nEquation, nTotalNode);

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            n_count[iNode] = 0;
            for (int m = 0; m < nEquation; ++ m)
            {
                q_n_list[m][iNode] = 0.0;
            }
        }

        typedef DataStruct_AdtTree<int, RDouble> AdtTree;
        AdtTree *adtTree;
        AdtTree::AdtNodeList nodeList;
        RDouble tol;
        CreateStandardADT(grid, adtTree, tol);
        RDouble coor[3], minwin[3], maxwin[3];

        RDouble3D &xs = *grid->GetStructX();
        RDouble3D &ys = *grid->GetStructY();
        RDouble3D &zs = *grid->GetStructZ();

        int count = 0;
        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            BCType = structBC->GetBCType();
            //if ((!(ISInterfaceBC(BCType))) && (!IsWall(true, BCType, 0))) continue;
            if (!(IsInterface(BCType)))
            {
                continue;
            }

            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            iSurface = structBC->GetFaceDirection() + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        structBC->GetInsideCellIndex(i, j, k, is, js, ks, 1);

                        RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);

                        structBC->GetGhostCellIndex(i, j, k, it, jt, kt, 1);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            qtmp[m] = half * (q(is, js, ks, m) + q(it, jt, kt, m));
                        }

                        structBC->GetBoundaryFaceIndex(i, j, k, il, jl, kl);

                        ir = irot[iSurface - 1];
                        jr = jrot[iSurface - 1];
                        kr = krot[iSurface - 1];

                        for (int m = 0; m < nlist; ++ m)
                        {
                            in = il + ijk[m][ir];
                            jn = jl + ijk[m][jr];
                            kn = kl + ijk[m][kr];

                            coor[0] = xs(in, jn, kn);
                            coor[1] = ys(in, jn, kn);
                            coor[2] = zs(in, jn, kn);

                            minwin[0] = coor[0] - tol;
                            minwin[1] = coor[1] - tol;
                            minwin[2] = coor[2] - tol;

                            maxwin[0] = coor[0] + tol;
                            maxwin[1] = coor[1] + tol;
                            maxwin[2] = coor[2] + tol;
                            GetCoorIndex(adtTree, coor, minwin, maxwin, count, index);

                            n_count[index] ++;

                            for (int mm = 0; mm < nEquation; ++ mm)
                            {
                                q_n_list[mm][index] += qtmp[mm];
                            }
                        }
                    }
                }
            }
        }

        int nINode = count;

        for (int i = 0; i < nINode; ++ i)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                q_n_list[m][i] /= n_count[i];
            }
        }

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            BCType = structBC->GetBCType();
            if (!(IsInterface(BCType)))
            {
                continue;
            }

            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            iSurface = structBC->GetFaceDirection() + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        structBC->GetBoundaryFaceIndex(i, j, k, il, jl, kl);

                        ir = irot[iSurface - 1];
                        jr = jrot[iSurface - 1];
                        kr = krot[iSurface - 1];

                        for (int m = 0; m < nlist; ++ m)
                        {
                            int in = il + ijk[m][ir];
                            int jn = jl + ijk[m][jr];
                            int kn = kl + ijk[m][kr];

                            coor[0] = xs(in, jn, kn);
                            coor[1] = ys(in, jn, kn);
                            coor[2] = zs(in, jn, kn);

                            minwin[0] = coor[0] - tol;
                            minwin[1] = coor[1] - tol;
                            minwin[2] = coor[2] - tol;

                            maxwin[0] = coor[0] + tol;
                            maxwin[1] = coor[1] + tol;
                            maxwin[2] = coor[2] + tol;

                            GetCoorIndex(adtTree, coor, minwin, maxwin, count, index);

                            for (int mm = 0; mm < nEquation; ++ mm)
                            {
                                q_n(in, jn, kn, mm) = q_n_list[mm][index];
                            }
                        }
                    }
                }
            }
        }

        delete adtTree;
        delete[]n_count;
        delete[]qtmp;
        DelPointer2(q_n_list);
    }

    void NSSolverStruct::CreateStandardADT(Grid *grid, DataStruct_AdtTree<int, RDouble> *&adtTree, RDouble &tol)
    {
        typedef DataStruct_AdtTree<int, RDouble> AdtTree;
        RDouble mindis, maxdis;
        grid->GetMinMaxDS(mindis, maxdis);

        tol = mindis / 4;

        RDouble pmin[3], pmax[3];

        RDouble *ptmin = grid->GetMinBox();
        RDouble *ptmax = grid->GetMaxBox();

        pmin[0] = ptmin[0];
        pmin[1] = ptmin[1];
        pmin[2] = ptmin[2];

        pmax[0] = ptmax[0];
        pmax[1] = ptmax[1];
        pmax[2] = ptmax[2];

        ShiftMinMaxBox(pmin, pmax, two * tol);

        adtTree = new AdtTree(3, pmin, pmax);
    }

    void NSSolverStruct::WallValueFix(Grid *gridIn, FieldProxy *qnProxy, FieldProxy *qProxy)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &q = qProxy->GetField_STR();
        RDouble4D &q_n = qnProxy->GetField_STR();

        int nEquation = GetNumberOfEquations();

        int nTotalNode = grid->GetNTotalNode();

        int nlist = 4;
        if (grid->GetDim() == TWO_D)
        {
            nlist = 2;
        }

        int ijk[4][3];
        int irot[3], jrot[3], krot[3];

        GetIJKROT(ijk, irot, jrot, krot);

        int ist, ied, jst, jed, kst, ked;
        int iSurface, BCType, index;
        int il, jl, kl;
        int is, js, ks;
        int it, jt, kt;
        int ir, jr, kr;

        RDouble *qtmp = new RDouble[nEquation];
        int *n_count = new int[nTotalNode];
        RDouble **q_n_list = NewPointer2<RDouble>(nEquation, nTotalNode);

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            n_count[iNode] = 0;
            for (int m = 0; m < nEquation; ++ m)
            {
                q_n_list[m][iNode] = 0.0;
            }
        }

        typedef DataStruct_AdtTree<int, RDouble> AdtTree;
        AdtTree *adtTree;
        AdtTree::AdtNodeList nodeList;
        RDouble tol;
        CreateStandardADT(grid, adtTree, tol);
        RDouble coor[3], minwin[3], maxwin[3];

        RDouble3D &xs = *grid->GetStructX();
        RDouble3D &ys = *grid->GetStructY();
        RDouble3D &zs = *grid->GetStructZ();

        int count = 0;
        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            BCType = structBC->GetBCType();
            if ((!IsWall(BCType)))
            {
                continue;
            }

            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            iSurface = structBC->GetFaceDirection() + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        structBC->GetInsideCellIndex(i, j, k, is, js, ks, 1);

                        RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);

                        structBC->GetGhostCellIndex(i, j, k, it, jt, kt, 1);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            qtmp[m] = half * (q(is, js, ks, m) + q(it, jt, kt, m));
                        }

                        structBC->GetBoundaryFaceIndex(i, j, k, il, jl, kl);

                        ir = irot[iSurface - 1];
                        jr = jrot[iSurface - 1];
                        kr = krot[iSurface - 1];

                        for (int m = 0; m < nlist; ++ m)
                        {
                            int in = il + ijk[m][ir];
                            int jn = jl + ijk[m][jr];
                            int kn = kl + ijk[m][kr];

                            coor[0] = xs(in, jn, kn);
                            coor[1] = ys(in, jn, kn);
                            coor[2] = zs(in, jn, kn);

                            minwin[0] = coor[0] - tol;
                            minwin[1] = coor[1] - tol;
                            minwin[2] = coor[2] - tol;

                            maxwin[0] = coor[0] + tol;
                            maxwin[1] = coor[1] + tol;
                            maxwin[2] = coor[2] + tol;
                            GetCoorIndex(adtTree, coor, minwin, maxwin, count, index);
                            n_count[index] ++;
                            for (int mm = 0; mm < nEquation; ++ mm)
                            {
                                q_n_list[mm][index] += qtmp[mm];
                            }
                        }
                    }
                }
            }
        }

        int nINode = count;

        for (int i = 0; i < nINode; ++ i)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                q_n_list[m][i] /= n_count[i];
            }
        }

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            BCType = structBC->GetBCType();
            if ((!IsWall(BCType)))
            {
                continue;
            }

            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            iSurface = structBC->GetFaceDirection() + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        structBC->GetBoundaryFaceIndex(i, j, k, il, jl, kl);

                        ir = irot[iSurface - 1];
                        jr = jrot[iSurface - 1];
                        kr = krot[iSurface - 1];

                        for (int m = 0; m < nlist; ++ m)
                        {
                            int in = il + ijk[m][ir];
                            int jn = jl + ijk[m][jr];
                            int kn = kl + ijk[m][kr];

                            coor[0] = xs(in, jn, kn);
                            coor[1] = ys(in, jn, kn);
                            coor[2] = zs(in, jn, kn);

                            minwin[0] = coor[0] - tol;
                            minwin[1] = coor[1] - tol;
                            minwin[2] = coor[2] - tol;

                            maxwin[0] = coor[0] + tol;
                            maxwin[1] = coor[1] + tol;
                            maxwin[2] = coor[2] + tol;

                            GetCoorIndex(adtTree, coor, minwin, maxwin, count, index);

                            for (int mm = 0; mm < nEquation; ++ mm)
                            {
                                q_n(in, jn, kn, mm) = q_n_list[mm][index];
                            }
                        }
                    }
                }
            }
        }

        delete adtTree;
        delete[]n_count;
        delete[]qtmp;
        DelPointer2(q_n_list);
    }

    void NSSolverStruct::InitFlowAsRestart()
    {
        StructGrid *grid = StructGridCast(GetGrid());

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        Param_NSSolverStruct *parameters = GetControlParameters();

        int outputStepInterval = 0;
        GlobalDataBase::UpdateData("outnstep", &outputStepInterval, PHINT, 1);

        RDouble physicalTime = 0.0;
        GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);

        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble4D &heatTransferCoeff = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("heatTransferCoeff"));

        RDouble refLenth = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

        RDouble3D &Xm = *(grid->GetCellCenterX());
        RDouble3D &Ym = *(grid->GetCellCenterY());
        RDouble3D &Zm = *(grid->GetCellCenterZ());

        int nEquation = GetNumberOfEquations();
        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

        RDouble iniSpeedCoef = parameters->GetIniSpeedCoef();

        Range M(0, nEquation - 1);
        int mst = M.first();
        int med = M.last();

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        primitiveVariables(i, j, k, m) = primitiveVarFarfield[m];
                    }
                }
            }
        }

        temperatures = 1.0;

#ifdef USE_HYFLOW_CODE
        using namespace PHENGLEI;
        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();

            if (BCType != INFLOW && BCType != PRESSURE_INLET)
            {
                continue;
            }

            RDouble *primitiveVarBC = new RDouble[nEquation];
            int initFlowBcType;
            RDouble firstPoint[3], secondPoint[3], bcRadius;

            //ini initFlowBcType
            if(!structBC->CheckParamData("initFlowBcType"))continue;
            structBC->GetParamData("primitiveVarFarfield", primitiveVarBC, PHDOUBLE, nEquation);
            structBC->GetParamData("initFlowBcType", &initFlowBcType, PHINT, 1);

            //ini firstPoint
            if(!structBC->CheckParamData("firstPoint"))continue;
            structBC->GetParamData("firstPoint",  &firstPoint , PHDOUBLE, 3);

            //ini secondPoint
            if(structBC->CheckParamData("secondPoint"))
            {
                structBC->GetParamData("secondPoint", &secondPoint, PHDOUBLE, 3);
            }
            else
            {
                secondPoint[0] = firstPoint[0];
                secondPoint[1] = firstPoint[1];
                secondPoint[2] = firstPoint[2];
            }
            
            //ini bcRadius
            if(structBC->CheckParamData("bcRadius"))
            {
                structBC->GetParamData("bcRadius", &bcRadius, PHDOUBLE, 1);
            }
            else
            {
                bcRadius = 0.0;
            }

            for (int s = 0; s < 3; ++ s)
            {
                firstPoint[s] /= refLenth;
                secondPoint[s] /= refLenth;
            }
            bcRadius /= refLenth;

            RDouble TBC[3];
            gas->GetTemperatureR(primitiveVarBC, TBC[0], TBC[1], TBC[2]);

            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {

                        if(isInside(Xm(i, j, k), Ym(i, j, k), Zm(i, j, k), firstPoint, secondPoint, bcRadius,initFlowBcType))
                        {
                            for (int m = mst; m <= med; ++ m)
                            {
                                primitiveVariables(i, j, k, m) = primitiveVarBC[m];
                            }

                            for (int m = 0; m < parameters->GetTemperatureModel(); ++ m)
                            {
                                temperatures(i, j, k, m) = TBC[m];
                            }
                        }
                    }
                }
            }
            delete[] primitiveVarBC;
        }

        if(iniSpeedCoef < 1.0 && iniSpeedCoef >= 0.0)
        {
            using namespace IDX;
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        primitiveVariables(i, j, k, IU) *= iniSpeedCoef;
                        primitiveVariables(i, j, k, IV) *= iniSpeedCoef;
                        primitiveVariables(i, j, k, IW) *= iniSpeedCoef;
                    }
                }
            }
        }
#endif

        InitGrad(grid);

        int isUnsteady = parameters->GetIsUnsteady();
        if (!isUnsteady)
        {
            RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
            residual = 0.0;
            return;
        }

        RDouble4D &primitiveVariablen1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
        RDouble4D &primitiveVariablen2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));

        primitiveVariablen1 = primitiveVariables;
        primitiveVariablen2 = primitiveVariablen1;

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &residualn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n1"));
        RDouble4D &residualn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n2"));

        residual = 0.0;
        residualn1 = 0.0;
        residualn2 = 0.0;

        int nStatisticalStep = 0;
        GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al.,
        //! exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble4D &resTmp = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_tmp"));

        resTmp = 0.0;
    }

    void NSSolverStruct::InitGrad(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &gradPrimtiveVarFaceX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
        RDouble4D &gradPrimtiveVarFaceY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
        RDouble4D &gradPrimtiveVarFaceZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));
        gradPrimtiveVarFaceX = 0.0;
        gradPrimtiveVarFaceY = 0.0;
        gradPrimtiveVarFaceZ = 0.0;

        RDouble4D &gradTemperatureFaceX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceX"));
        RDouble4D &gradTemperatureFaceY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceY"));
        RDouble4D &gradTemperatureFaceZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceZ"));
        gradTemperatureFaceX = 0.0;
        gradTemperatureFaceY = 0.0;
        gradTemperatureFaceZ = 0.0;

        RDouble4D &gradUVWTCellCenterX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
        RDouble4D &gradUVWTCellCenterY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
        RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));
        gradUVWTCellCenterX = 0.0;
        gradUVWTCellCenterY = 0.0;
        gradUVWTCellCenterZ = 0.0;
    }

    void NSSolverStruct::InitCGrid(Grid *fineGrid, Grid *coarseGrid)
    {

    }

    void NSSolverStruct::Boundary(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        Param_NSSolverStruct *parameters = GetControlParameters();
        int nChemical = parameters->GetChemicalFlag();
        int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

        int viscousType;
        RDouble wallTemperature;

        using namespace PHENGLEI;

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        int iWallBC = 0;
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();

            if (BCType < 0)
            {
                continue;
            }
            else if (BCType == EXTRAPOLATION)
            {
                OutFlowBC3D(grid, structBC);
            }
            else if (IsWall(BCType))
            {
                structBC->GetParamData("viscousType", &viscousType, PHINT, 1);
                structBC->GetParamData("wallTemperature", &wallTemperature, PHDOUBLE, 1);

                if (viscousType == INVISCID)
                {
                    SymmetryBC3D(grid, structBC);
                }
                else
                {
#ifdef USE_HYFLOW_CODE
                    //! Compute the temperature of wall.
                    GAS_SPACE::ComputeSurfaceTemperature(grid, structBC, iWallBC);

                    //! Compute the mass fractions of species on wall.
                    GAS_SPACE::ComputeSurfaceMassFraction(grid, structBC, iWallBC);

                    //! Compute the slip variables on wall.
                    GAS_SPACE::ComputeSurfaceSlipVariables(grid, structBC, iWallBC);
#endif

                    if (wallTemperature < 0.0)
                    {
                        if (nChemical > 0)
                        {
#ifdef USE_HYFLOW_CODE
                            GAS_SPACE::ChemicalVisWallISO3D(grid, structBC, iWallBC);
#endif
                        }
                        else    //! Viscous adiabatic wall condition.
                        {
                            ViscousAdiabaticWall(grid, structBC);
                        }
                    }
                    else//! Viscous iso-thermal wall condition.
                    {
                        if (nChemical > 0)
                        {
#ifdef USE_HYFLOW_CODE
                            GAS_SPACE::ChemicalVisWallISO3D(grid, structBC, iWallBC);
                            //GAS_SPACE::SetChemicalIsothermalWall(grid, StructBC);
#endif
                        }
                        else
                        {
                            ViscousIsotropicWall(grid, structBC, iWallBC);
                        }
                    }
                }
                ++ iWallBC;
            }
            else if (BCType == ABLATION_SURFACE)
            {
#ifdef USE_HYFLOW_CODE
                structBC->GetParamData("viscousType", &viscousType, PHINT, 1);
                if (viscousType > INVISCID)
                {
                    //! Compute the temperature and mass fractions of species under the ablation wall condition.(added by meijie in 2021.09)
                    GAS_SPACE::ComputeSurfaceAblation(grid, structBC, iWallBC);
                }

                GAS_SPACE::ChemicalAblationWall(grid, structBC, iWallBC);
#endif
                ++ iWallBC;
            }
            else if (BCType == SYMMETRY)
            {
                SymmetryBC3D(grid, structBC);
            }
            else if (BCType == FARFIELD)
            {
                if (ifLowSpeedPrecon == 0)
                {
                    FarFieldRiemannInvariants(grid, structBC);
                }
                else
                {
                    PreconFarfieldBC(grid, structBC);
                }
            }
            else if (BCType == INFLOW)
            {
                InFlowBC3D(grid, structBC);
            }
            else if (BCType == OUTFLOW)
            {
                OutFlowBC3D(grid, structBC);
            }
            else if (BCType == POLE || BCType / 10 == POLE)
            {
                SingularAxisBC3D(grid, structBC);
            }
            else if (BCType == GENERIC_1)
            {
                //! Do not deal with this condition now!
                continue;
            }
            else if (BCType == PRESSURE_INLET)
            {
                PressureInletBCRiemann(grid, structBC);
            }
            else if (BCType == PRESSURE_OUTLET)
            {
                PressureOutletBC(grid, structBC);
            }
            else if (BCType == MASS_FLOW_INLET)
            {
                MassFlowInletBC(grid, structBC);
            }
            else if (BCType == MASS_FLOW_OUTLET)
            {
                MassFlowOutletBC(grid, structBC);
            }
            else if (BCType == EXTERNAL_BC)
            {
                OutFlowBC3D(grid, structBC);
            }
            else if (BCType == NOZZLE_INFLOW)
            {
                NozzleInletFlow(grid, structBC);
            }
            else if (BCType == NOZZLE_OUTFLOW)
            {
                NozzleOutFlow(grid, structBC);
            }
            else if (BCType == OVERSET)
            {
                OutFlowBC3D(grid, structBC);
            }
            else
            {
                //! In order to avoid illegal BCtype ID.
                ostringstream oss;
                oss << "Error : Illegal BCtype ID " << BCType << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }
        //! The corner point is only a temporary solution and should be replaced as much as possible.
        CornerPoint(grid);
    }

    void NSSolverStruct::CornerPoint(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int nEquation = GetNumberOfEquations();

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        FillCornerPoint3D(q, ni, nj, nk, nEquation);
    }

    void NSSolverStruct::calculateBCArea(Grid *gridIn, StructBC *bcRegion)
    {
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        //BC total area does not change with calculating, so it only calculates once.
        if (bcData->IsExist("BCTotalArea", PHDOUBLE, 1))
        {
            return;
        }

        using namespace PHMPI;
        string bcName = bcRegion->GetBCName();
        RDouble BCTotalArea = 0.0;

        //one processor may owe some zones, so calculate all BC area in one processor.
        int nLocalZones = GetNumberofLocalZones();
        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int iZoneID = GetLocalZoneIDToGlobalZoneID(iZone);
            Grid *iGrid = PHSPACE::GetGrid(iZoneID);

            StructGrid *grid = StructGridCast(iGrid);

            CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
            int nBCRegion = compositeBCRegion->GetnBCRegion();
            for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
            {
                StructBC *jBCRegion = compositeBCRegion->GetBCRegion(iBCRegion);
                string jBCName = jBCRegion->GetBCName();

                //judge the same BC region according the BC name.
                if (jBCName == bcName)
                {
                    int ist, ied, jst, jed, kst, ked;
                    int in, jn, kn;
                    jBCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);
                    int iSurface = jBCRegion->GetFaceDirection() + 1;
                    RDouble4D &area = *(grid->GetFaceArea());

                    for (int k = kst; k <= ked; ++ k)
                    {
                        for (int j = jst; j <= jed; ++ j)
                        {
                            for (int i = ist; i <= ied; ++ i)
                            {
                                jBCRegion->GetBoundaryFaceIndex(i, j, k, in, jn, kn);
                                BCTotalArea += area(in, jn, kn, iSurface);
                            }
                        }
                    }
                }
            }
        }

        //get the total area in all processors with the same BC name.
        RDouble GLB_BCTotalArea;
        PH_AllReduce(&BCTotalArea, &GLB_BCTotalArea, 1, MPI_SUM);
        BCTotalArea = GLB_BCTotalArea;

        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int iZoneID = GetLocalZoneIDToGlobalZoneID(iZone);
            Grid *iGrid = PHSPACE::GetGrid(iZoneID);

            StructGrid *grid = StructGridCast(iGrid);

            CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
            int nBCRegion = compositeBCRegion->GetnBCRegion();
            for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
            {
                StructBC *jBCRegion = compositeBCRegion->GetBCRegion(iBCRegion);
                string jBCName = jBCRegion->GetBCName();

                if (jBCName == bcName)
                {
                    Data_Param *bcData = jBCRegion->GetBCParamDataBase();
                    bcData->UpdateData("BCTotalArea", &BCTotalArea, PHDOUBLE, 1);
                }
            }
        }
    }

    void NSSolverStruct::calculateMassFluxRatio(Grid *gridIn, StructBC *bcRegion)
    {
        using namespace PHMPI;
        string BCName = bcRegion->GetBCName();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();
        RDouble massFlow = 0.0;
        bcData->GetData("massFlow", &massFlow, PHDOUBLE, 1);
        RDouble massFlux = 0.0;

        Param_NSSolverStruct *parameters = GetControlParameters();
        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

        //one processor may owe some zones, so calculate all BC mass flux in one processor.
        //in fact, massflux need calculte once in one step, but I dont find a good method to judge if it is
        //calculated or not, so it will calculate more times, may cost a little more time.
        int nLocalZones = GetNumberofLocalZones();
        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int iZoneID = GetLocalZoneIDToGlobalZoneID(iZone);
            Grid *iGrid = PHSPACE::GetGrid(iZoneID);
            StructGrid *grid = StructGridCast(iGrid);

            CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
            int nBCRegion = compositeBCRegion->GetnBCRegion();
            for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
            {
                StructBC *jBCRegion = compositeBCRegion->GetBCRegion(iBCRegion);
                string jBCName = jBCRegion->GetBCName();

                //judge the same BC region according the BC name.
                if (jBCName == BCName)
                {
                    int ist, ied, jst, jed, kst, ked;
                    int in, jn, kn;
                    int is, js, ks;
                    RDouble nxs, nys, nzs;
                    int layer = 1;
                    jBCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);
                    int iSurface = jBCRegion->GetFaceDirection() + 1;
                    int lr = jBCRegion->GetFaceLeftOrRightIndex();
                    RDouble4D &area = *(grid->GetFaceArea());
                    RDouble4D &xfn = *(grid->GetFaceNormalX());
                    RDouble4D &yfn = *(grid->GetFaceNormalY());
                    RDouble4D &zfn = *(grid->GetFaceNormalZ());
                    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
                    using namespace IDX;
                    for (int k = kst; k <= ked; ++ k)
                    {
                        for (int j = jst; j <= jed; ++ j)
                        {
                            for (int i = ist; i <= ied; ++ i)
                            {
                                jBCRegion->GetInsideCellIndex(i, j, k, is, js, ks, layer);
                                jBCRegion->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                                if (ABS(q(is, js, ks, IU)) > LARGE)
                                {
                                    q(is, js, ks, IR) = primitiveVarFarfield[IR];
                                    q(is, js, ks, IU) = primitiveVarFarfield[IU];
                                    q(is, js, ks, IV) = primitiveVarFarfield[IV];
                                    q(is, js, ks, IW) = primitiveVarFarfield[IW];
                                }

                                nxs = lr * xfn(in, jn, kn, iSurface);
                                nys = lr * yfn(in, jn, kn, iSurface);
                                nzs = lr * zfn(in, jn, kn, iSurface);
                                RDouble vn_c1 = q(is, js, ks, IU) * nxs + q(is, js, ks, IV) * nys + q(is, js, ks, IW) * nzs;

                                massFlux += q(is, js, ks, IR) * vn_c1 * area(in, jn, kn, iSurface);
                            }
                        }
                    }
                }
            }
        }

        //! Get the total area in all processors with the same BC name.
        RDouble GLB_MassFlux;
        PH_AllReduce(&massFlux, &GLB_MassFlux, 1, MPI_SUM);

        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
        RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();

        //calculte the dimensional mass flux in this BC name region.
        GLB_MassFlux *= refDimensionalDensity * refDimensionalSonicSpeed * refMachNumber;
        RDouble massFluxRatio = massFlow / (GLB_MassFlux + TINY);

        massFluxRatio = MAX(0.8, MIN(massFluxRatio, 1.2));
        bcData->UpdateData("massFluxRatio", &massFluxRatio, PHDOUBLE, 1);
    }

    void NSSolverStruct::OutFlowBC3D(Grid *gridIn, StructBC *structBC)
    {
        using namespace IDX;
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int is1, js1, ks1, is2, js2, ks2;
        int it1, jt1, kt1, it2, jt2, kt2;

        Param_NSSolverStruct *parameters = GetControlParameters();
        int isWennScheme = parameters->GetWennSchemeFlag();

        int nEquation = GetNumberOfEquations();

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    //! First cell index inside flowfield.
                    structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);

                    //! Second cell index inside flowfield.
                    structBC->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                    RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                    RestrictIndex(is2, js2, ks2, ni - 1, nj - 1, nk - 1);

                    //! First cell index of ghostcell.
                    structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);

                    //! Second cell index of ghostcell.
                    structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        q(it1, jt1, kt1, m) = q(is1, js1, ks1, m);
                        q(it2, jt2, kt2, m) = q(it1, jt1, kt1, m);
                    }

                    if (isWennScheme == 1)
                    {
                        int is3, js3, ks3;
                        int it3, jt3, kt3;

                        //! Third cell index inside flow field.
                        structBC->GetInsideCellIndex(i, j, k, is3, js3, ks3, 3);

                        RestrictIndex(is3, js3, ks3, ni - 1, nj - 1, nk - 1);

                        //! Third cell index of ghost cell.
                        structBC->GetGhostCellIndex(i, j, k, it3, jt3, kt3, 3);

                        for (int m = 0; m < nEquation; ++m)
                        {
                            q(it3, jt3, kt3, m) = q(it1, jt1, kt1, m);
                        }
                    }
                }
            }
        }
    }

    void NSSolverStruct::SymmetryBC3D(Grid *gridIn, StructBC *structBC)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &faceNormalVelocity = *(grid->GetFaceNormalVelocity());

        int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        int iSurface = structBC->GetFaceDirection() + 1;

        RDouble nxs = 0, nys = 0, nzs = 0;

        int is = 0, js = 0, ks = 0, it = 0, jt = 0, kt = 0, in = 0, jn = 0, kn = 0;
        RDouble vx = 0, vy = 0, vz = 0, vn = 0, vgn = 0;
        int nLayer = GetNumberOfGhostCellLayers();

        using namespace IDX;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                    nxs = xfn(in, jn, kn, iSurface);
                    nys = yfn(in, jn, kn, iSurface);
                    nzs = zfn(in, jn, kn, iSurface);
                    vgn = faceNormalVelocity(in, jn, kn, iSurface);

                    for (int layer = 1; layer <= nLayer; ++ layer)
                    {
                        structBC->GetInsideCellIndex(i, j, k, is, js, ks, layer);
                        RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);
                        structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                        vx = q(is, js, ks, IU);
                        vy = q(is, js, ks, IV);
                        vz = q(is, js, ks, IW);
                        vn = nxs * vx + nys * vy + nzs * vz - vgn;

                        q(it, jt, kt, IR) = q(is, js, ks, IR);
                        q(it, jt, kt, IU) = q(is, js, ks, IU) - 2.0 * nxs * vn;
                        q(it, jt, kt, IV) = q(is, js, ks, IV) - 2.0 * nys * vn;
                        q(it, jt, kt, IW) = q(is, js, ks, IW) - 2.0 * nzs * vn;
                        q(it, jt, kt, IP) = q(is, js, ks, IP);

                        for (int m = nNSEquation; m < nEquation; ++ m)
                        {
                            q(it, jt, kt, m) = q(is, js, ks, m);
                        }
                    }
                }
            }
        }
    }

    void NSSolverStruct::SingularAxisBC3D(Grid *gridIn, StructBC *structBC)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &faceNormalVelocity = *(grid->GetFaceNormalVelocity());

        int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        int iSurface = structBC->GetFaceDirection() + 1;

        int iii[3][3] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
        //! To judge if the singular axis using half-plane symmetry boundary conditions or using full plane boundary conditions
        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        int iSingularDim = 0;
        bool isWhole = 1;
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();
            if (BCType == PHENGLEI::SYMMETRY)
            {
                isWhole = 0;
                iSingularDim = structBC->GetFaceDirection() + 1;
                if (iSingularDim != iSurface) break;
            }

            if (BCType < 0)
            {
                if (structBC->GetRegionBlock() == structBC->GetTargetRegionBlock())
                {
                    if (structBC->GetFaceDirectionOfTargetBlock() == iSingularDim)
                    {
                        isWhole = 1;
                        iSingularDim = structBC->GetFaceDirection() + 1;
                        if (iSingularDim != iSurface) break;
                    }
                }
            }
        }

        int *startPoint = structBC->GetStartPoint();
        int *endPoint = structBC->GetEndPoint();
        int *faceDirectionIndex = structBC->GetFaceDirectionIndex();
        if (iSingularDim == 0)
        {
            SingularAxisFullPlaneBC3D(grid, structBC);
            return;
        }
        int nd3 = 6 - iSingularDim - iSurface;
        RDouble nxs = 0, nys = 0, nzs = 0, vgn = 0;
        if (isWhole == 0)
        {
            if (iSingularDim == 1)
            {
                int i = startPoint[0];
                int j = (endPoint[1] + startPoint[1]) / 2 - 3 * faceDirectionIndex[1];
                int k = (endPoint[2] + startPoint[2]) / 2 - 3 * faceDirectionIndex[2];
                nxs = xfn(i, j, k, iSingularDim);
                nys = yfn(i, j, k, iSingularDim);
                nzs = zfn(i, j, k, iSingularDim);
                vgn = faceNormalVelocity(i, j, k, iSingularDim);
            }
            else if (iSingularDim == 2)
            {
                int i = (endPoint[0] + startPoint[0]) / 2 - 3 * faceDirectionIndex[0];
                int j = startPoint[1];
                int k = (endPoint[2] + startPoint[2]) / 2 - 3 * faceDirectionIndex[2];
                nxs = xfn(i, j, k, iSingularDim);
                nys = yfn(i, j, k, iSingularDim);
                nzs = zfn(i, j, k, iSingularDim);
                vgn = faceNormalVelocity(i, j, k, iSingularDim);
            }
            else if (iSingularDim == 3)
            {
                int i = (endPoint[0] + startPoint[0]) / 2 - 3 * faceDirectionIndex[0];
                int j = (endPoint[1] + startPoint[1]) / 2 - 3 * faceDirectionIndex[1];
                int k = startPoint[2];
                nxs = xfn(i, j, k, iSingularDim);
                nys = yfn(i, j, k, iSingularDim);
                nzs = zfn(i, j, k, iSingularDim);
                vgn = faceNormalVelocity(i, j, k, iSingularDim);
            }
        }
        else
        {
            nxs = 0.0;
            nys = 0.0;
            nzs = 0.0;
            vgn = 0.0;
        }

        int halfDim = startPoint[iSingularDim - 1] + (endPoint[iSingularDim - 1] - startPoint[iSingularDim - 1]) / (1 + isWhole);
        int wholeDim = endPoint[iSingularDim - 1];

        RDouble vx = 0, vy = 0, vz = 0, vn = 0;

        using namespace IDX;

        for (int i = startPoint[iSurface - 1]; i <= endPoint[iSurface - 1]; ++ i)
        {
            for (int iLayer = 1; iLayer <= 2; ++ iLayer)
            {
                int it1 = i + iLayer * structBC->GetFaceLeftOrRightIndex();
                int is1 = i - (iLayer - 1) * structBC->GetFaceLeftOrRightIndex();
                for (int j = startPoint[iSingularDim - 1]; j <= endPoint[iSingularDim - 1]; ++ j)
                {
                    int jj = 0;
                    if (isWhole == 0)
                    {
                        jj = wholeDim - j + 1;
                    }
                    else
                    {
                        jj = halfDim + j;
                        if (jj > wholeDim) jj = wholeDim;
                    }

                    for (int k = startPoint[nd3 - 1]; k <= endPoint[nd3 - 1]; ++ k)
                    {
                        int it = iii[iSurface - 1][0] * it1 + iii[iSingularDim - 1][0] * j + iii[nd3 - 1][0] * k;
                        int jt = iii[iSurface - 1][1] * it1 + iii[iSingularDim - 1][1] * j + iii[nd3 - 1][1] * k;
                        int kt = iii[iSurface - 1][2] * it1 + iii[iSingularDim - 1][2] * j + iii[nd3 - 1][2] * k;
                        int is = iii[iSurface - 1][0] * is1 + iii[iSingularDim - 1][0] * jj + iii[nd3 - 1][0] * k;
                        int js = iii[iSurface - 1][1] * is1 + iii[iSingularDim - 1][1] * jj + iii[nd3 - 1][1] * k;
                        int ks = iii[iSurface - 1][2] * is1 + iii[iSingularDim - 1][2] * jj + iii[nd3 - 1][2] * k;

                        vx = q(is, js, ks, IU);
                        vy = q(is, js, ks, IV);
                        vz = q(is, js, ks, IW);
                        vn = nxs * vx + nys * vy + nzs * vz - vgn;

                        q(it, jt, kt, IR) = q(is, js, ks, IR);
                        q(it, jt, kt, IU) = q(is, js, ks, IU) - 2.0 * nxs * vn;
                        q(it, jt, kt, IV) = q(is, js, ks, IV) - 2.0 * nys * vn;
                        q(it, jt, kt, IW) = q(is, js, ks, IW) - 2.0 * nzs * vn;
                        q(it, jt, kt, IP) = q(is, js, ks, IP);

                        for (int m = nNSEquation; m < nEquation; ++ m)
                        {
                            q(it, jt, kt, m) = q(is, js, ks, m);
                        }
                    }
                }
            }
        }
    }

    void NSSolverStruct::SingularAxisFullPlaneBC3D(Grid *gridIn, StructBC *structBC)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        int bcType = structBC->GetBCType();
        int SingularDim = bcType - 71;

        int is1, js1, ks1;
        int it1, jt1, kt1, it2, jt2, kt2;

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        if (SingularDim == 0)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int iEq = 0; iEq < nEquation; ++ iEq)
                {
                    RDouble sumQ = 0.0;
                    if (jst == jed)
                    {
                        int j = jst;
                        for (int k = kst; k <= ked; ++ k)
                        {
                            structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                            RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                            sumQ += q(is1, js1, ks1, iEq);
                        }
                        sumQ = sumQ / static_cast<RDouble>(ked - kst + 1);
                        for (int k = kst; k <= ked; ++ k)
                        {
                            structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                            structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);
                            q(it1, jt1, kt1, iEq) = sumQ;
                            q(it2, jt2, kt2, iEq) = q(it1, jt1, kt1, iEq);
                        }
                    }
                    else
                    {
                        int k = kst;
                        for (int j = jst; j <= jed; ++ j)
                        {
                            structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                            RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                            sumQ += q(is1, js1, ks1, iEq);
                        }
                        sumQ = sumQ / static_cast<RDouble>(jed - jst + 1);
                        for (int j = jst; j <= ked; ++ j)
                        {
                            structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                            structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);
                            q(it1, jt1, kt1, iEq) = sumQ;
                            q(it2, jt2, kt2, iEq) = q(it1, jt1, kt1, iEq);
                        }
                    }
                }
            }
        }
        else if (SingularDim == 1)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int iEq = 0; iEq < nEquation; ++ iEq)
                {
                    RDouble sumQ = 0.0;
                    if (ist == ied)
                    {
                        int i = ist;
                        for (int k = kst; k <= ked; ++ k)
                        {
                            structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                            sumQ += q(is1, js1, ks1, iEq);
                        }
                        sumQ = sumQ / static_cast<RDouble>(ked - kst + 1);
                        for (int k = kst; k <= ked; ++ k)
                        {
                            structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                            structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);
                            q(it1, jt1, kt1, iEq) = sumQ;
                            q(it2, jt2, kt2, iEq) = q(it1, jt1, kt1, iEq);
                        }
                    }
                    else
                    {
                        int k = kst;
                        for (int i = ist; i <= ied; ++ i)
                        {
                            structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                            sumQ += q(is1, js1, ks1, iEq);
                        }
                        sumQ = sumQ / static_cast<RDouble>(ied - ist + 1);
                        for (int i = ist; i <= ied; ++ i)
                        {
                            structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                            structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);
                            q(it1, jt1, kt1, iEq) = sumQ;
                            q(it2, jt2, kt2, iEq) = q(it1, jt1, kt1, iEq);
                        }
                    }
                }
            }
        }
        else if (SingularDim == 2)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int iEq = 0; iEq < nEquation; ++ iEq)
                {
                    RDouble sumQ = 0.0;
                    if (jst == jed)
                    {
                        int j = jst;
                        for (int i = ist; i <= ied; ++ i)
                        {
                            structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                            RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                            sumQ += q(is1, js1, ks1, iEq);
                        }
                        sumQ = sumQ / static_cast<RDouble>(ied - ist + 1);
                        for (int i = ist; i <= ied; ++ i)
                        {
                            structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                            structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);
                            q(it1, jt1, kt1, iEq) = sumQ;
                            q(it2, jt2, kt2, iEq) = q(it1, jt1, kt1, iEq);
                        }
                    }
                    else
                    {
                        int i = ist;
                        for (int j = jst; j <= jed; ++ j)
                        {
                            structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                            RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                            sumQ += q(is1, js1, ks1, iEq);
                        }
                        sumQ = sumQ / static_cast<RDouble>(jed - jst + 1);
                        for (int j = jst; j <= ked; ++ j)
                        {
                            structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                            structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);
                            q(it1, jt1, kt1, iEq) = sumQ;
                            q(it2, jt2, kt2, iEq) = q(it1, jt1, kt1, iEq);
                        }
                    }
                }
            }
        }
    }

    //! Set viscous adiabatic wall condition.
    void NSSolverStruct::ViscousAdiabaticWall(Grid *gridIn, StructBC *structBC)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int isWennScheme = parameters->GetWennSchemeFlag();
        //! Obtain the number of the physical equations.
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

        RDouble4D &faceVelocityX = *(grid->GetFaceVelocityX());
        RDouble4D &faceVelocityY = *(grid->GetFaceVelocityY());
        RDouble4D &faceVelocityZ = *(grid->GetFaceVelocityZ());

        RDouble uWall = 0.0;
        RDouble vWall = 0.0;
        RDouble wWall = 0.0;
        Data_Param *bcData = structBC->GetBCParamDataBase();

        if (bcData)
        {
            if (bcData->IsExist("uWall", PHDOUBLE, 1))
            {
                bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
                bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
                bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
            }
        }

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        int iSurface = structBC->GetFaceDirection() + 1;

        int is1, js1, ks1, is2, js2, ks2;
        int it1, jt1, kt1, it2, jt2, kt2;

        using namespace IDX;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int in, jn, kn;
                    structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                    structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                    structBC->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                    RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                    RestrictIndex(is2, js2, ks2, ni - 1, nj - 1, nk - 1);

                    structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                    structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                    RDouble velocityXWall = uWall;
                    RDouble velocityYWall = vWall;
                    RDouble velocityZWall = wWall;
                    if (isUnsteady && isAle)
                    {
                        velocityXWall += faceVelocityX(in, jn, kn, iSurface);
                        velocityYWall += faceVelocityY(in, jn, kn, iSurface);
                        velocityZWall += faceVelocityZ(in, jn, kn, iSurface);
                    }

                    q(it1, jt1, kt1, IR) = q(is1, js1, ks1, IR);
                    q(it1, jt1, kt1, IU) = -q(is1, js1, ks1, IU) + 2.0 * velocityXWall;
                    q(it1, jt1, kt1, IV) = -q(is1, js1, ks1, IV) + 2.0 * velocityYWall;
                    q(it1, jt1, kt1, IW) = -q(is1, js1, ks1, IW) + 2.0 * velocityZWall;
                    q(it1, jt1, kt1, IP) = q(is1, js1, ks1, IP);

                    q(it2, jt2, kt2, IR) = q(it1, jt1, kt1, IR);
                    q(it2, jt2, kt2, IU) = q(it1, jt1, kt1, IU);
                    q(it2, jt2, kt2, IV) = q(it1, jt1, kt1, IV);
                    q(it2, jt2, kt2, IW) = q(it1, jt1, kt1, IW);
                    q(it2, jt2, kt2, IP) = q(it1, jt1, kt1, IP);

                    for (int m = nNSEquation; m < nEquation; ++ m)
                    {
                        q(it1, jt1, kt1, m) = q(is1, js1, ks1, m);
                        q(it2, jt2, kt2, m) = q(it1, jt1, kt1, m);
                    }

                    if (isWennScheme == 1)
                    {
                        int is3, js3, ks3;
                        int it3, jt3, kt3;

                        //! Third cell index inside flowfield.
                        structBC->GetInsideCellIndex(i, j, k, is3, js3, ks3, 3);

                        RestrictIndex(is3, js3, ks3, ni - 1, nj - 1, nk - 1);

                        //! Third cell index of ghostcell.
                        structBC->GetGhostCellIndex(i, j, k, it3, jt3, kt3, 3);

                        q(it3, jt3, kt3, IR) = q(it1, jt1, kt1, IR);
                        q(it3, jt3, kt3, IU) = q(it1, jt1, kt1, IU);
                        q(it3, jt3, kt3, IV) = q(it1, jt1, kt1, IV);
                        q(it3, jt3, kt3, IW) = q(it1, jt1, kt1, IW);
                        q(it3, jt3, kt3, IP) = q(it1, jt1, kt1, IP);

                        for (int m = 0; m < nEquation; ++m)
                        {
                            q(it3, jt3, kt3, m) = q(it1, jt1, kt1, m);
                        }
                    }
                }
            }
        }
    }

    //! Set viscous isothermal wall condition.
    void NSSolverStruct::ViscousIsotropicWall(Grid *gridIn, StructBC *structBC, int iWallBC)
    {
        //! chemical = 1 is not need to be cared, it is treated in another function.
        //! temperature array should not appear here, it has not been updated.
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Param_NSSolverStruct *parameters = GetControlParameters();
        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

        RDouble refGama = parameters->GetRefGama();
        RDouble refMachNumber = parameters->GetRefMachNumber();

        RDouble uWall = 0.0;
        RDouble vWall = 0.0;
        RDouble wWall = 0.0;

        RDouble wallTemperature = 0.0;
         structBC->GetParamData("wallTemperature", &wallTemperature, PHDOUBLE, 1);

        if (structBC->CheckParamData("uWall"))
        {
            structBC->GetParamData("uWall", &uWall, PHDOUBLE, 1);
            structBC->GetParamData("vWall", &vWall, PHDOUBLE, 1);
            structBC->GetParamData("wWall", &wWall, PHDOUBLE, 1);
        }

        bool wallTempArrayExist = false;
        RDouble3D *wallTempArray = NULL;
        if (structBC->CheckFieldData("wallTempArray"))
        {
            wallTempArrayExist = true;
            wallTempArray = reinterpret_cast <RDouble3D *> (structBC->GetFieldDataPtr("wallTempArray"));
        }
        RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        RDouble3D *surfaceSlipVariables = NULL;
        if (nSlipBCModel > 0)
        {
            surfaceSlipVariables = reinterpret_cast <RDouble3D *> (gridIn->GetDataPtr("surfaceSlipVariables"));
        }
        RDouble2D *surfaceTemperature = NULL;
        if (wallTemperature >= 0.0)
        {
            surfaceTemperature = reinterpret_cast <RDouble2D *> (gridIn->GetDataPtr("surfaceTemperature"));
        }
        bool isRadiationEquilibrium = false;
        if (fabs(wallTemperature) <= EPSILON)
        {
            isRadiationEquilibrium = true;
        }

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);
        int iSurface = structBC->GetFaceDirection() + 1;

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &faceVelocityX = *(grid->GetFaceVelocityX());
        RDouble4D &faceVelocityY = *(grid->GetFaceVelocityY());
        RDouble4D &faceVelocityZ = *(grid->GetFaceVelocityZ());

        RDouble temperatureWallNonDimensional = wallTemperature / refDimensionalTemperature;  //! non-dimensional temperature.
        RDouble temperaturelimitation = 1.E-6;

        int nLayer = GetNumberOfGhostCellLayers();

        using namespace IDX;

        int is, js, ks, it, jt, kt;
        int nIndexOfCell = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int layer = 1; layer <= nLayer; ++ layer)
                    {
                        int in, jn, kn;
                        structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                        structBC->GetInsideCellIndex(i, j, k, is, js, ks, layer);
                        RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);
                        structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                        RDouble rs = q(is, js, ks, 0);
                        RDouble us = q(is, js, ks, 1);
                        RDouble vs = q(is, js, ks, 2);
                        RDouble ws = q(is, js, ks, 3);
                        RDouble ps = q(is, js, ks, 4);
                        RDouble ts = refGama * refMachNumber * refMachNumber * ps / rs;

                        if (nSlipBCModel > 0)
                        {
                            uWall = (*surfaceSlipVariables)(iWallBC, nIndexOfCell, 3);
                            vWall = (*surfaceSlipVariables)(iWallBC, nIndexOfCell, 4);
                            wWall = (*surfaceSlipVariables)(iWallBC, nIndexOfCell, 5);
                        }

                        RDouble velocityXWall = uWall;
                        RDouble velocityYWall = vWall;
                        RDouble velocityZWall = wWall;
                        if (isUnsteady && isAle)
                        {
                            velocityXWall += faceVelocityX(in, jn, kn, iSurface);
                            velocityYWall += faceVelocityY(in, jn, kn, iSurface);
                            velocityZWall += faceVelocityZ(in, jn, kn, iSurface);
                        }

                        RDouble ut = -us + 2.0 * velocityXWall;
                        RDouble vt = -vs + 2.0 * velocityYWall;
                        RDouble wt = -ws + 2.0 * velocityZWall;
                        RDouble pt = ps;
                        if (wallTempArrayExist)
                        {
                            temperatureWallNonDimensional = (*wallTempArray)(i, j, k) / refDimensionalTemperature;
                        }
                        if (isRadiationEquilibrium)
                        {
                            temperatureWallNonDimensional = (*surfaceTemperature)(iWallBC, nIndexOfCell);
                        }
                        if (nSlipBCModel > 0)
                        {
                            temperatureWallNonDimensional = (*surfaceSlipVariables)(iWallBC, nIndexOfCell, ITT);
                        }

                        RDouble tt = 2.0 * temperatureWallNonDimensional - ts;
                        if (tt < temperaturelimitation)
                        {
                            tt = temperatureWallNonDimensional;
                        }

                        RDouble rt = refGama * refMachNumber * refMachNumber * pt / tt;

                        q(it, jt, kt, 0) = rt;
                        q(it, jt, kt, 1) = ut;
                        q(it, jt, kt, 2) = vt;
                        q(it, jt, kt, 3) = wt;
                        q(it, jt, kt, 4) = pt;
                    }
                    ++ nIndexOfCell;
                }
            }
        }
    }

    void NSSolverStruct::SetBCInitialValuesByInterpolation(StructBC *structBC, RDouble *primitiveVar)
    {
        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        if (!isUnsteady)
        {
            return;
        }

        using namespace GAS_SPACE;
        //Obtain the local parameters on boundary.
        int dataSize = structBC->GetBCParamDataBase()->GetSizeOfData("time");
        RDouble *timeVector = new RDouble[dataSize];
        structBC->GetBCParamDataBase()->GetData("time", &timeVector[0], PHDOUBLE, dataSize);
        RDouble *pressureVector = new RDouble[dataSize];
        structBC->GetBCParamDataBase()->GetData("pressure", &timeVector[0], PHDOUBLE, dataSize);
        RDouble *velocityVector = new RDouble[dataSize];
        structBC->GetBCParamDataBase()->GetData("velocity", &timeVector[0], PHDOUBLE, dataSize);
        RDouble *temperatureVector = new RDouble[dataSize];
        structBC->GetBCParamDataBase()->GetData("temperature", &timeVector[0], PHDOUBLE, dataSize);

        RDouble localDimensionalVelocity = 0, localDimensionalTemperature = 0, localDimensionalPressure = 0, localDimensionalDensity = 0;

        //Obtain the reference parameters.
        RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
        RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

        RDouble physicalTime = GlobalDataBase::GetDoubleParaFromDB("physicalTime") / refDimensionalVelocity;
        RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
        RDouble unit = reynoldsReferenceLengthDimensional / refDimensionalVelocity;
        RDouble physicalTimeStepDimensional = physicalTime * unit;

        localDimensionalVelocity    =    velocityVector[0];
        localDimensionalPressure    =    pressureVector[0];
        localDimensionalTemperature = temperatureVector[0];
        for (int iData = 1; iData < dataSize; ++ iData)
        {
            if (physicalTimeStepDimensional > timeVector[iData - 1] && physicalTimeStepDimensional < timeVector[iData])
            {
                RDouble delta = (timeVector[iData] - physicalTimeStepDimensional) / (timeVector[iData] - timeVector[iData - 1]);
                localDimensionalVelocity    =    velocityVector[iData - 1] + delta * (   velocityVector[iData] -    velocityVector[iData - 1]);
                localDimensionalPressure    =    pressureVector[iData - 1] + delta * (   pressureVector[iData] -    pressureVector[iData - 1]);
                localDimensionalTemperature = temperatureVector[iData - 1] + delta * (temperatureVector[iData] - temperatureVector[iData - 1]);
            }
        }
        if(physicalTimeStepDimensional > timeVector[dataSize - 1])
        {
            localDimensionalVelocity    =    velocityVector[dataSize - 1];
            localDimensionalPressure    =    pressureVector[dataSize - 1];
            localDimensionalTemperature = temperatureVector[dataSize - 1];
        }

        RDouble attackd = 0, angleSlide = 0;
        if (structBC->GetBCParamDataBase()->IsExist("attackd", PHDOUBLE, 1))
        {
            structBC->GetBCParamDataBase()->GetData("attackd", &attackd, PHDOUBLE, 1);
        }
        else
        {
            attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
        }
        if (structBC->GetBCParamDataBase()->IsExist("angleSlide", PHDOUBLE, 1))
        {
            structBC->GetBCParamDataBase()->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
        }
        RDouble attack   = attackd * PI / 180.0;
        RDouble sideslip = angleSlide * PI / 180.0;

        int nchem = gas->GetChemicalType();
        if (nchem == 0)
        {
            RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");
            localDimensionalDensity = localDimensionalPressure / reference_general_gas_constant / localDimensionalTemperature;
        }
        else
        {
            int nNSEqn = gas->GetNSEquationNumber();
            int numberOfSpecies = gas->GetNumberOfSpecies();
            RDouble initMassFraction[MAX_SPECIES_NUM] = { 0 }, primVars[MAX_SPECIES_NUM] = { 0 };
            if (structBC->GetBCParamDataBase()->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
            {
                structBC->GetBCParamDataBase()->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            }
            else
            {
                RDouble *speciesMass = gas->GetInitMassFraction();
                for (int s = 0; s < numberOfSpecies; ++ s)
                {
                    initMassFraction[s] = speciesMass[s];
                }
            }
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                primVars[nNSEqn + s] = initMassFraction[s];
            }

            RDouble oMass = 0.0, gasR = 0.0;
            gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
            gasR = rjmk * oMass;

            localDimensionalDensity = localDimensionalPressure / (gasR * localDimensionalTemperature);

            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                primitiveVar[nNSEqn + s] = initMassFraction[s];
            }
        }

        using namespace IDX;
        primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;
        primitiveVar[IU] = localDimensionalVelocity / refDimensionalVelocity * cos(attack) * cos(sideslip);
        primitiveVar[IV] = localDimensionalVelocity / refDimensionalVelocity * sin(attack) * cos(sideslip);
        primitiveVar[IW] = localDimensionalVelocity / refDimensionalVelocity * sin(sideslip);
        primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);
        delete []timeVector; timeVector = NULL;
        delete []pressureVector; pressureVector = NULL;
        delete []velocityVector; velocityVector = NULL;
        delete []temperatureVector; temperatureVector = NULL;
    }

    void NSSolverStruct::InFlowBC3D(Grid *gridIn, StructBC *structBC)
    {
        StructGrid *grid = StructGridCast(gridIn);

        Param_NSSolverStruct *parameters = GetControlParameters();

        int nEquation = GetNumberOfEquations();

        RDouble4D &primitiveVars = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        RDouble *primitiveVarFarfield = new RDouble[nEquation];
        structBC->GetParamData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

        int inflowParaType = -1;
        if (structBC->CheckParamData("inflowParaType"))
        {
            structBC->GetParamData("inflowParaType", &inflowParaType, PHINT, 1);
        }

        if (inflowParaType == TRAJECTORY)
        {
            SetBCInitialValuesByInterpolation(structBC, primitiveVarFarfield);
        }

        int it, jt, kt;

        int nLayer = GetNumberOfGhostCellLayers();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int layer = 1; layer <= nLayer; ++ layer)
                    {
                        structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            primitiveVars(it, jt, kt, m) = primitiveVarFarfield[m];
                        }
                    }
                }
            }
        }
        delete[] primitiveVarFarfield;
    }

    //! He Xin. Numerical simulation and qualitative analysis investigation of the supersonic flows over sharp cones at high angle of attack[D].
    //! China Aerodynamics Research and Development Center graduate department, December 2002.
    void NSSolverStruct::FarFieldRiemannInvariants(Grid *gridIn, StructBC *structBC)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int isWennScheme = parameters->GetWennSchemeFlag();

        int nEquation = GetNumberOfEquations();

        RDouble *primitiveVarFarfield = new RDouble[nEquation];
        structBC->GetParamData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());
        RDouble4D &faceNormalVelocity = *(grid->GetFaceNormalVelocity());

        RDouble refGama = parameters->GetRefGama();

        //! s_lr3d is used for calculating id, jd, kd, and the latters are further used for calculating in, jn, kn.
        //! in, jn, kn, iSurface and lr are used for calculating nx, ny, nz.
        int iSurface = structBC->GetFaceDirection() + 1;
        int lr = structBC->GetFaceLeftOrRightIndex();

        int in, jn, kn;
        RDouble nxs, nys, nzs;

        //! The index of infinite condition.
        using namespace IDX;
        RDouble roo, uoo, voo, woo, poo, vnoo, coo; //! vnoo is the projection on (nx,ny,nz).
        roo = primitiveVarFarfield[IR];
        uoo = primitiveVarFarfield[IU];
        voo = primitiveVarFarfield[IV];
        woo = primitiveVarFarfield[IW];
        poo = primitiveVarFarfield[IP];

        //! The index of inner, ghost and boundary points.
        int is1, js1, ks1, is2, js2, ks2;    //! The first and second inner points.
        int it1, jt1, kt1, it2, jt2, kt2;    //! The first and second ghost points.

        RDouble *prims1, *prims2;    //! The inner point.
        RDouble *primt1, *primt2;    //! The ghost point.
        prims1 = new RDouble[nEquation];
        prims2 = new RDouble[nEquation];
        primt1 = new RDouble[nEquation];
        primt2 = new RDouble[nEquation];

        RDouble rin, uin, vin, win, pin, vnin, cinner, vein; //! inner point variables.
        RDouble rb, ub, vb, wb, pb, vnb, cb;           //! boundary point variables.

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                    structBC->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                    //! The function of "RestrictIndex" has no use generally, and it may enforce the wrong value to the unknown value and lead to unexpected calculation.
                    //! However, it need more tests to decide whether delete it here or not.
                    RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                    RestrictIndex(is2, js2, ks2, ni - 1, nj - 1, nk - 1);

                    structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                    structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                    //! Initialization.
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prims1[m] = q(is1, js1, ks1, m);
                        prims2[m] = q(is2, js2, ks2, m);

                        primt1[m] = prims1[m];
                        primt2[m] = prims2[m];
                    }

                    //! The needed (nx, ny, nz) is stored in the cell of (in, jn, kn).
                    structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                    nxs = lr * xfn(in, jn, kn, iSurface);
                    nys = lr * yfn(in, jn, kn, iSurface);
                    nzs = lr * zfn(in, jn, kn, iSurface);
                    RDouble vgn = faceNormalVelocity(in, jn, kn, iSurface);

                    //! The first inner point.
                    rin = prims1[IR];
                    uin = prims1[IU];
                    vin = prims1[IV];
                    win = prims1[IW];
                    pin = prims1[IP];
                    vnin = nxs * uin + nys * vin + nzs * win - vgn;
                    cinner = sqrt(ABS(refGama * pin / rin));
                    vein = sqrt(uin * uin + vin * vin + win * win);

                    //! Infinite.
                    vnoo = nxs * uoo + nys * voo + nzs * woo - vgn;
                    coo = sqrt(ABS(refGama * poo / roo));

                    if (vein > cinner)
                    {
                        //! Supersonic.
                        if (vnin >= 0.0)
                        {
                            //! Supersonic outlet.
                            for (int m = 0; m < nEquation; ++ m)
                            {
                                q(it1, jt1, kt1, m) = prims1[m];
                                q(it2, jt2, kt2, m) = 2.0 * prims1[m] - prims2[m];
                            }

                            if (q(it2, jt2, kt2, IR) <= 0.0 || q(it2, jt2, kt2, IP) <= 0.0)
                            {
                                for (int m = 0; m < nEquation; ++ m)
                                {
                                    q(it2, jt2, kt2, m) = prims1[m];
                                }
                            }
                        }
                        else
                        {
                            //! Supersonic inlet.
                            for (int m = 0; m < nEquation; ++ m)
                            {
                                q(it1, jt1, kt1, m) = primitiveVarFarfield[m];
                                q(it2, jt2, kt2, m) = primitiveVarFarfield[m];
                            }
                        }
                    }
                    else
                    {
                        //! Subsonic.
                        RDouble gama1 = refGama - 1.0;
                        RDouble riemp;     //! Riemann invariant for inner.
                        RDouble riemm;     //! Riemann invariant for infinite.
                        riemp = vnin + 2.0 * cinner / gama1;
                        riemm = vnoo - 2.0 * coo / gama1;

                        vnb = half * (riemp + riemm);
                        cb = fourth * (riemp - riemm) * gama1;

                        RDouble sb;                          //! Entropy at the boundary.
                        RDouble uref, vref, wref, vnref;     //! Outlet refers to inner point, and inlet refers to infinite.

                        if (vnb > 0.0)
                        {
                            //! Subsonic outlet.
                            sb = pin / pow(rin, refGama);
                            uref = uin;
                            vref = vin;
                            wref = win;
                            vnref = vnin;
                        }
                        else
                        {
                            //! Subsonic inlet.
                            sb = poo / pow(roo, refGama);
                            uref = uoo;
                            vref = voo;
                            wref = woo;
                            vnref = vnoo;
                        }

                        rb = pow((cb * cb / (sb * refGama)), 1.0 / gama1);
                        ub = uref + nxs * (vnb - vnref);
                        vb = vref + nys * (vnb - vnref);
                        wb = wref + nzs * (vnb - vnref);
                        pb = cb * cb * rb / refGama;     //! The same to "pb = sb * pow(rb, refGama);".

                        primt1[IR] = rb;
                        primt1[IU] = ub;
                        primt1[IV] = vb;
                        primt1[IW] = wb;
                        primt1[IP] = pb;

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            primt2[m] = primt1[m];
                        }

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            q(it1, jt1, kt1, m) = primt1[m];
                            q(it2, jt2, kt2, m) = primt2[m];
                        }
                    }
                }
            }
        }

        if (isWennScheme == 1)
        {
            int is3, js3, ks3;
            int it3, jt3, kt3;
            RDouble* prims3;
            RDouble* primt3;
            prims3 = new RDouble[nEquation];
            primt3 = new RDouble[nEquation];

            for (int k = kst; k <= ked; ++k)
            {
                for (int j = jst; j <= jed; ++j)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        //! Third cell index inside flowfield.
                        structBC->GetInsideCellIndex(i, j, k, is3, js3, ks3, 3);

                        RestrictIndex(is3, js3, ks3, ni - 1, nj - 1, nk - 1);

                        //! Third cell index of ghostcell.
                        structBC->GetGhostCellIndex(i, j, k, it3, jt3, kt3, 3);

                        for (int m = 0; m < nEquation; ++m)
                        {
                            prims3[m] = q(is3, js3, ks3, m);
                            primt3[m] = prims3[m];
                        }
                        //! The needed (nx, ny, nz) is stored in the cell of (in, jn, kn).
                        structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                        nxs = lr * xfn(in, jn, kn, iSurface);
                        nys = lr * yfn(in, jn, kn, iSurface);
                        nzs = lr * zfn(in, jn, kn, iSurface);
                        RDouble vgn = faceNormalVelocity(in, jn, kn, iSurface);

                        //! The first inner point.
                        rin = prims1[IR];
                        uin = prims1[IU];
                        vin = prims1[IV];
                        win = prims1[IW];
                        pin = prims1[IP];
                        vnin = nxs * uin + nys * vin + nzs * win - vgn;
                        cinner = sqrt(ABS(refGama * pin / rin));
                        vein = sqrt(uin * uin + vin * vin + win * win);

                        //! Infinite.
                        vnoo = nxs * uoo + nys * voo + nzs * woo - vgn;
                        coo = sqrt(ABS(refGama * poo / roo));

                        if (vein > cinner)
                        {
                            //! Supersonic.
                            if (vnin >= 0.0)
                            {
                                //! Supersonic outlet.
                                for (int m = 0; m < nEquation; ++m)
                                {
                                    q(it3, jt3, kt3, m) = 3.0 * prims1[m] - 2.0 * prims2[m];
                                }

                                if (q(it3, jt3, kt3, IR) <= 0.0 || q(it3, jt3, kt3, IP) <= 0.0)
                                {
                                    for (int m = 0; m < nEquation; ++m)
                                    {
                                        q(it3, jt3, kt3, m) = prims1[m];
                                    }
                                }
                            }
                            else
                            {
                                //! Supersonic inlet.
                                for (int m = 0; m < nEquation; ++m)
                                {
                                    q(it3, jt3, kt3, m) = primitiveVarFarfield[m];
                                }
                            }
                        }
                        else
                        {
                            //! Subsonic.
                            RDouble gama1 = refGama - 1.0;
                            RDouble riemp;     //! Riemann invariant for inner.
                            RDouble riemm;     //! Riemann invariant for infinite.
                            riemp = vnin + 2.0 * cinner / gama1;
                            riemm = vnoo - 2.0 * coo / gama1;

                            vnb = half * (riemp + riemm);
                            cb = fourth * (riemp - riemm) * gama1;

                            RDouble sb;                          //! Entropy at the boundary.
                            RDouble uref, vref, wref, vnref;     //! Outlet refers to inner point, and inlet refers to infinite.

                            if (vnb > 0.0)
                            {
                                //! Subsonic outlet.
                                sb = pin / pow(rin, refGama);
                                uref = uin;
                                vref = vin;
                                wref = win;
                                vnref = vnin;
                            }
                            else
                            {
                                //! Subsonic inlet.
                                sb = poo / pow(roo, refGama);
                                uref = uoo;
                                vref = voo;
                                wref = woo;
                                vnref = vnoo;
                            }

                            rb = pow((cb * cb / (sb * refGama)), 1.0 / gama1);
                            ub = uref + nxs * (vnb - vnref);
                            vb = vref + nys * (vnb - vnref);
                            wb = wref + nzs * (vnb - vnref);
                            pb = cb * cb * rb / refGama;     //! The same to "pb = sb * pow(rb, refGama);".

                            primt1[IR] = rb;
                            primt1[IU] = ub;
                            primt1[IV] = vb;
                            primt1[IW] = wb;
                            primt1[IP] = pb;

                            for (int m = 0; m < nEquation; ++m)
                            {
                                primt3[m] = primt1[m];
                            }

                            for (int m = 0; m < nEquation; ++m)
                            {
                                q(it3, jt3, kt3, m) = primt3[m];
                            }
                        }
                    }
                }
            }

            delete[] prims3;
            delete[] primt3;
        }

        delete[] prims1;
        delete[] prims2;
        delete[] primt1;
        delete[] primt2;

        delete[] primitiveVarFarfield;
    }

    void NSSolverStruct::PreconFarfieldBC(Grid *gridIn, StructBC *structBC)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int isWennScheme = parameters->GetWennSchemeFlag();

        int nEquation = GetNumberOfEquations();

        RDouble *primitiveVarFarfield = new RDouble[nEquation];
        structBC->GetParamData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());

        RDouble refGama = parameters->GetRefGama();
        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble mach2 = refMachNumber * refMachNumber;
        //! in, jn, kn, iSurface and lr are used for calculating nx, ny, nz.
        int iSurface = structBC->GetFaceDirection() + 1;
        int lr = structBC->GetFaceLeftOrRightIndex();

        int in, jn, kn;
        RDouble nxs, nys, nzs;

        //! The index of infinite condition.
        using namespace IDX;
        RDouble roo, uoo, voo, woo, poo, too, vnoo, coo; //! vnoo is the projection on (nx,ny,nz).
        roo = primitiveVarFarfield[IR];
        uoo = primitiveVarFarfield[IU];
        voo = primitiveVarFarfield[IV];
        woo = primitiveVarFarfield[IW];
        poo = primitiveVarFarfield[IP];
        too = refGama * mach2 * poo / roo;
        //! The index of inner, ghost and boundary points.
        int is1, js1, ks1, is2, js2, ks2;   //! The first and second inner points.
        int it1, jt1, kt1, it2, jt2, kt2;   //! The first and second ghost points.

        RDouble *prims1, *prims2;  //! The inner point.
        RDouble *primt1, *primt2;  //! The ghost point.
        prims1 = new RDouble[nEquation];
        prims2 = new RDouble[nEquation];
        primt1 = new RDouble[nEquation];
        primt2 = new RDouble[nEquation];

        RDouble rin, uin, vin, win, pin, tin, vnin, cinner, vein;    //! Inner point variables.
        RDouble rb, ub, vb, wb, pb, tb;    //! Boundary point variables.

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                    structBC->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                    //! The function of "RestrictIndex" has no use generally, and it may enforce the wrong value to the unknown value and lead to unexpected calculation.
                    //! However, it need more tests to decide whether delete it here or not.
                    RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                    RestrictIndex(is2, js2, ks2, ni - 1, nj - 1, nk - 1);

                    structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                    structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                    //! Initialization.
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prims1[m] = q(is1, js1, ks1, m);
                        prims2[m] = q(is2, js2, ks2, m);

                        primt1[m] = prims1[m];
                        primt2[m] = prims2[m];
                    }

                    //! The needed (nx, ny, nz) is stored in the cell of (in, jn, kn).
                    structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                    nxs = lr * xfn(in, jn, kn, iSurface);
                    nys = lr * yfn(in, jn, kn, iSurface);
                    nzs = lr * zfn(in, jn, kn, iSurface);

                    //! The first inner point.
                    rin = prims1[IR];
                    uin = prims1[IU];
                    vin = prims1[IV];
                    win = prims1[IW];
                    pin = prims1[IP];
                    tin = refGama * mach2 * pin / rin;
                    vnin = nxs * uin + nys * vin + nzs * win; //! Not consider the "vgn".
                    cinner = sqrt(ABS(refGama * pin / rin));
                    vein = sqrt(uin * uin + vin * vin + win * win);
                    RDouble rmach = ABS(vein / cinner);
                    //! Infinite.
                    vnoo = nxs * uoo + nys * voo + nzs * woo;
                    coo = sqrt(ABS(refGama * poo / roo));

                    if (rmach < 1.0)
                    {
                        //! Subsonic.
                        if (vnin < 0.0)
                        {
                            ub = uoo;
                            vb = voo;
                            wb = woo;
                            pb = pin;
                            tb = too;
                            rb = refGama * mach2 * pb / tb;
                        }
                        else
                        {
                            ub = uin;
                            vb = vin;
                            wb = win;
                            pb = poo;
                            tb = tin;
                            rb = refGama * mach2 * pb / tb;
                        }

                        primt1[IR] = rb;
                        primt1[IU] = ub;
                        primt1[IV] = vb;
                        primt1[IW] = wb;
                        primt1[IP] = pb;

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            primt2[m] = primt1[m];
                        }

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            q(it1, jt1, kt1, m) = primt1[m];
                            q(it2, jt2, kt2, m) = primt2[m];
                        }
                    }
                    else
                    {
                        //! Supersonic.
                        if (vnin < 0.0)
                        {
                            //! Supersonic inlet.
                            for (int m = 0; m < nEquation; ++ m)
                            {
                                q(it1, jt1, kt1, m) = primitiveVarFarfield[m];
                                q(it2, jt2, kt2, m) = primitiveVarFarfield[m];
                            }
                        }
                        else
                        {
                            //! Supersonic outlet.
                            for (int m = 0; m < nEquation; ++ m)
                            {
                                q(it1, jt1, kt1, m) = prims1[m];
                                q(it2, jt2, kt2, m) = prims1[m];
                            }
                        }
                    }
                }
            }
        }

        if (isWennScheme == 1)
        {
            int is3, js3, ks3;
            int it3, jt3, kt3;
            RDouble* prims3;
            RDouble* primt3;
            prims3 = new RDouble[nEquation];
            primt3 = new RDouble[nEquation];

            for (int k = kst; k <= ked; ++k)
            {
                for (int j = jst; j <= jed; ++j)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        //! Third cell index inside flowfield.
                        structBC->GetInsideCellIndex(i, j, k, is3, js3, ks3, 3);

                        RestrictIndex(is3, js3, ks3, ni - 1, nj - 1, nk - 1);

                        //! Third cell index of ghostcell.
                        structBC->GetGhostCellIndex(i, j, k, it3, jt3, kt3, 3);

                        for (int m = 0; m < nEquation; ++m)
                        {
                            prims3[m] = q(is3, js3, ks3, m);
                            primt3[m] = prims3[m];
                        }

                        //! The needed (nx, ny, nz) is stored in the cell of (in, jn, kn).
                        structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                        nxs = lr * xfn(in, jn, kn, iSurface);
                        nys = lr * yfn(in, jn, kn, iSurface);
                        nzs = lr * zfn(in, jn, kn, iSurface);

                        //! The first inner point.
                        rin = prims1[IR];
                        uin = prims1[IU];
                        vin = prims1[IV];
                        win = prims1[IW];
                        pin = prims1[IP];
                        tin = refGama * mach2 * pin / rin;
                        vnin = nxs * uin + nys * vin + nzs * win; //! Not consider the "vgn".
                        cinner = sqrt(ABS(refGama * pin / rin));
                        vein = sqrt(uin * uin + vin * vin + win * win);
                        RDouble rmach = ABS(vein / cinner);
                        //! Infinite.
                        vnoo = nxs * uoo + nys * voo + nzs * woo;
                        coo = sqrt(ABS(refGama * poo / roo));

                        if (rmach < 1.0)
                        {
                            //! Subsonic.
                            if (vnin < 0.0)
                            {
                                ub = uoo;
                                vb = voo;
                                wb = woo;
                                pb = pin;
                                tb = too;
                                rb = refGama * mach2 * pb / tb;
                            }
                            else
                            {
                                ub = uin;
                                vb = vin;
                                wb = win;
                                pb = poo;
                                tb = tin;
                                rb = refGama * mach2 * pb / tb;
                            }

                            primt1[IR] = rb;
                            primt1[IU] = ub;
                            primt1[IV] = vb;
                            primt1[IW] = wb;
                            primt1[IP] = pb;

                            for (int m = 0; m < nEquation; ++m)
                            {
                                primt3[m] = primt1[m];
                            }

                            for (int m = 0; m < nEquation; ++m)
                            {
                                q(it3, jt3, kt3, m) = primt3[m];
                            }
                        }
                        else
                        {
                            //! Supersonic.
                            if (vnin < 0.0)
                            {
                                //! Supersonic inlet.
                                for (int m = 0; m < nEquation; ++m)
                                {
                                    q(it3, jt3, kt3, m) = primitiveVarFarfield[m];
                                }
                            }
                            else
                            {
                                //! Supersonic outlet.
                                for (int m = 0; m < nEquation; ++m)
                                {
                                    q(it3, jt3, kt3, m) = prims1[m];
                                }
                            }
                        }
                    }
                }
            }
            delete[] prims3;
            delete[] primt3;
        }
        delete[] prims1;
        delete[] prims2;
        delete[] primt1;
        delete[] primt2;

        delete[] primitiveVarFarfield;
    }

    void NSSolverStruct::PressureInletBCRiemann(Grid *gridIn, StructBC *structBC)
    {
        using namespace IDX;

        Param_NSSolverStruct *parameters = GetControlParameters();
        int isWennScheme = parameters->GetWennSchemeFlag();

        int nEquation = GetNumberOfEquations();

        RDouble *primitiveVarFarfield = new RDouble[nEquation];
        Data_Param *bcData = structBC->GetBCParamDataBase();
        bcData->GetData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

        RDouble totalPressure, totalTemperature, direction_inlet[3];
        bcData->GetData("totalPressure", &totalPressure, PHDOUBLE, 1);
        bcData->GetData("totalTemperature", &totalTemperature, PHDOUBLE, 1);
        bcData->GetData("direction_inlet", &direction_inlet, PHDOUBLE, 3);

        //! Dimensional.
        RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
        RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
        RDouble velocity = refMachNumber * refDimensionalSonicSpeed;
        RDouble dynamic_pressure = refDimensionalDensity * velocity * velocity;
        RDouble refGama = parameters->GetRefGama();
        RDouble gama1 = refGama - 1;

        int nChemical = parameters->GetChemicalFlag();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        Thermo_Energy &Thermo_Energy_temparay = *gas->GetThermo_Energy_temparay();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nTemperatureModel = parameters->GetTemperatureModel();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0};

        int mTT = parameters->GetmTT();
        int mTV = parameters->GetmTV();
        int mTE = parameters->GetmTE();
        int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");

        if(nChemical > 0)
        {
            if (bcData->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
            {
                bcData->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
                for (int s = 0; s < numberOfSpecies; ++ s)
                {
                    primitiveVarFarfield[nNSEquation + s] = initMassFraction[s];
                }
            }
            else
            {
                RDouble *speciesMass = gas->GetInitMassFraction();
                for (int s = 0; s < numberOfSpecies; ++ s)
                {
                    initMassFraction[s] = speciesMass[s];
                }
            }
        }

        //! Nondimensionalization.
        totalPressure /= dynamic_pressure;
        totalTemperature /= refDimensionalTemperature;

        RDouble coefficientOfStateEquation = 1.0 / (refGama * refMachNumber * refMachNumber);

        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        int iSurface = structBC->GetFaceDirection() + 1;
        int lr = structBC->GetFaceLeftOrRightIndex();
        int is, js, ks, it, jt, kt;
        int in, jn, kn;
        RDouble nxs, nys, nzs;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    structBC->GetInsideCellIndex(i, j, k, is, js, ks, 1);
                    RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);

                    structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                    nxs = lr * xfn(in, jn, kn, iSurface);
                    nys = lr * yfn(in, jn, kn, iSurface);
                    nzs = lr * zfn(in, jn, kn, iSurface);

                    if(nChemical > 0)  //  for chemical  by dms
                    {
                        structBC->GetGhostCellIndex(i, j, k, it, jt, kt, 1);  //ghost
                        if(outIterStep > 1)
                        {
                            gas->ComputeMixturegasThermalParameter(initMassFraction, t(it, jt, kt,mTT), t(it, jt, kt, mTT), t(it, jt, kt, mTT), 0);   //  equlibrim
                        }
                        else
                        {
                            gas->ComputeMixturegasThermalParameter(initMassFraction, totalTemperature, totalTemperature, totalTemperature, 0);   //  equlibrim
                        }
                        
                        refGama = Thermo_Energy_temparay.gama;
                        gama1 = refGama - 1;
                        coefficientOfStateEquation = 1.0 / (refGama * refMachNumber * refMachNumber);

                        if(nTemperatureModel == 2)
                        {
                            primitiveVarFarfield[nEquation - 1] = Thermo_Energy_temparay.Ev + Thermo_Energy_temparay.Ee;
                        }
                        else if(nTemperatureModel == 3)
                        {
                            primitiveVarFarfield[nEquation - 1] = Thermo_Energy_temparay.Ee;
                            primitiveVarFarfield[nEquation - 2] = Thermo_Energy_temparay.Ev;
                        }
                    }

                    for (int layer = 1; layer <= 1; ++ layer)
                    {
                        structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                        //get inner cell flow field information.
                        RDouble rInside = q(is, js, ks, IR);
                        RDouble V2Inside = q(is, js, ks, IU) * q(is, js, ks, IU) + q(is, js, ks, IV) * q(is, js, ks, IV)
                            + q(is, js, ks, IW) * q(is, js, ks, IW);
                        RDouble pInside = q(is, js, ks, IP);

                        //calculate the total energy, total enthalpy and sound speed.
                        RDouble totalEnergy = pInside / (rInside * gama1) + 0.5 * V2Inside;
                        RDouble totalEnthalpy = (refGama * coefficientOfStateEquation / gama1) * totalTemperature;
                        RDouble c2 = refGama * pInside / rInside;

                        //calculate the riemann invariants R+.
                        RDouble Riemann = 2.0 * sqrt(c2) / gama1;
                        Riemann += q(is, js, ks, IU) * nxs + q(is, js, ks, IV) * nys + q(is, js, ks, IW) * nzs;

                        RDouble theta = nxs * direction_inlet[0] + nys * direction_inlet[1] + nzs * direction_inlet[2];

                        RDouble totalc2 = gama1 * (totalEnthalpy - (totalEnergy + pInside / rInside) + 0.5 * V2Inside) + c2;     //Co^2

                        RDouble a = 1.0 + 0.5 * gama1 * theta * theta;
                        RDouble b = -1.0 * gama1 * theta * Riemann;
                        RDouble c = 0.5 * gama1 * Riemann * Riemann - 2.0 * totalc2 / gama1;

                        RDouble d = b * b - 4.0 * a * c;
                        d = sqrt(max(0.0, d));
                        RDouble VInside = (-b + d) / (2.0 * a);
                        VInside = max(0.0, VInside);
                        V2Inside = VInside * VInside;

                        c2 = totalc2 - 0.5 * gama1 * V2Inside;

                        RDouble Mach2 = V2Inside / c2;
                        Mach2 = min(1.0, Mach2);
                        V2Inside = Mach2 * c2;
                        VInside = sqrt(V2Inside);
                        c2 = totalc2 - 0.5 * gama1 * V2Inside;

                        RDouble tInside = c2 / (refGama * coefficientOfStateEquation);
                        pInside = totalPressure * pow((tInside / totalTemperature), refGama / gama1);
                        rInside = pInside / (coefficientOfStateEquation * tInside);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            q(it, jt, kt, m) = primitiveVarFarfield[m];
                        }

                        q(it, jt, kt, IR) = rInside;
                        q(it, jt, kt, IU) = VInside * direction_inlet[0];
                        q(it, jt, kt, IV) = VInside * direction_inlet[1];
                        q(it, jt, kt, IW) = VInside * direction_inlet[2];
                        q(it, jt, kt, IP) = pInside;

                        for (int m = 0; m < nTemperatureModel; ++ m)
                        {
                            t(it, jt, kt, m) = tInside;
                        }
                    }

                    int it2, jt2, kt2;
                    structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        q(it2, jt2, kt2, m) = q(it, jt, kt, m);
                    }

                    if (isWennScheme == 1)
                    {
                        int it3, jt3, kt3;
                        structBC->GetGhostCellIndex(i, j, k, it3, jt3, kt3, 3);

                        for (int m = 0; m < nEquation; ++m)
                        {
                            q(it3, jt3, kt3, m) = q(it, jt, kt, m);
                        }
                    }
                }
            }
        }
        delete[] primitiveVarFarfield;
    }

    void NSSolverStruct::MassFlowInletBC(Grid *gridIn, StructBC *structBC)
    {
        calculateBCArea(gridIn, structBC);

        // get the basic mesh information.
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());

        // get the flow field information.
        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        int nEquation = GetNumberOfEquations();

        using namespace IDX;
        Data_Param *bcData = structBC->GetBCParamDataBase();

        // get the basic variables for mass flow in boundary.
        RDouble totalTemperature, massFlow, BCTotalArea, direction_inlet[3];
        bcData->GetData("massFlow", &massFlow, PHDOUBLE, 1);
        bcData->GetData("totalTemperature", &totalTemperature, PHDOUBLE, 1);
        bcData->GetData("direction_inlet", &direction_inlet, PHDOUBLE, 3);
        bcData->GetData("BCTotalArea", &BCTotalArea, PHDOUBLE, 1);

        //! Dimensional.
        Param_NSSolverStruct *parameters = GetControlParameters();
        RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
        RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
        RDouble refGama = parameters->GetRefGama();
        int isWennScheme = parameters->GetWennSchemeFlag();

        //! Nondimensionalization.
        massFlow /= BCTotalArea;
        totalTemperature /= refDimensionalTemperature;
        RDouble coefficientOfStateEquation = 1.0 / (refGama * refMachNumber * refMachNumber);

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        int iSurface = structBC->GetFaceDirection() + 1;
        int lr = structBC->GetFaceLeftOrRightIndex();
        int is, js, ks, it, jt, kt;
        int in, jn, kn;
        RDouble nxs, nys, nzs;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    structBC->GetInsideCellIndex(i, j, k, is, js, ks, 1);
                    RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);

                    structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                    nxs = lr * xfn(in, jn, kn, iSurface);
                    nys = lr * yfn(in, jn, kn, iSurface);
                    nzs = lr * zfn(in, jn, kn, iSurface);

                    for (int layer = 1; layer <= 1; ++ layer)
                    {
                        structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                        q(it, jt, kt, IR) = q(is, js, ks, IR);
                        RDouble Vnb = massFlow / (refDimensionalDensity * refDimensionalSonicSpeed * refMachNumber) / q(is, js, ks, IR);  //注意：此处massFlow为单位面积流量

                        q(it, jt, kt, IU) = Vnb * direction_inlet[0];
                        q(it, jt, kt, IV) = Vnb * direction_inlet[1];
                        q(it, jt, kt, IW) = Vnb * direction_inlet[2];

                        //square of the stagnation sound speed.
                        RDouble c2 = refGama * coefficientOfStateEquation * totalTemperature;
                        //square of the critical sound speed.
                        RDouble c2_ct = 2.0 * c2 / (refGama + 1.0);
                        RDouble lamd2 = Vnb * Vnb / c2_ct;
                        RDouble tb = totalTemperature * (1.0 - (refGama - 1.0) / (refGama + 1.0) * lamd2);
                        q(it, jt, kt, IP) = q(is, js, ks, IR) * coefficientOfStateEquation * tb;
                    }

                    int it2, jt2, kt2;
                    structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        q(it2, jt2, kt2, m) = q(it, jt, kt, m);
                    }

                    if (isWennScheme == 1)
                    {
                        int it3, jt3, kt3;
                        structBC->GetGhostCellIndex(i, j, k, it3, jt3, kt3, 3);
                        for (int m = 0; m < nEquation; ++m)
                        {
                            q(it3, jt3, kt3, m) = q(it, jt, kt, m);
                        }
                    }
                }
            }
        }
    }

    void NSSolverStruct::MassFlowOutletBC(Grid *gridIn, StructBC *structBC)
    {
        // calculate BC area and mass flux ratio of this BC region.
        calculateBCArea(gridIn, structBC);
        calculateMassFluxRatio(gridIn, structBC);

        using namespace IDX;
        //! Dimensional.
        Param_NSSolverStruct *parameters = GetControlParameters();
        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
        RDouble refGama = parameters->GetRefGama();

        Data_Param *bcData = structBC->GetBCParamDataBase();
        RDouble massFluxRatio, massFlow, BCTotalArea;
        bcData->GetData("massFluxRatio", &massFluxRatio, PHDOUBLE, 1);
        bcData->GetData("massFlow", &massFlow, PHDOUBLE, 1);
        bcData->GetData("BCTotalArea", &BCTotalArea, PHDOUBLE, 1);

        //! Nondimensionalization.
        RDouble coefficientOfStateEquation = 1.0 / (refGama * refMachNumber * refMachNumber);

        // fluxDensity : rho * u;
        RDouble fluxDensity = massFlow / (refDimensionalDensity * refDimensionalDensity * refMachNumber) / BCTotalArea;

        RDouble mfTbCoeff = 0.5 * (refGama - 1.0) / refGama / coefficientOfStateEquation *
            (1.0 - massFluxRatio * massFluxRatio);

        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        int iSurface = structBC->GetFaceDirection() + 1;
        int lr = structBC->GetFaceLeftOrRightIndex();
        int is, js, ks, it, jt, kt;
        int in, jn, kn;
        RDouble nxs, nys, nzs;

        int nLayer = GetNumberOfGhostCellLayers();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    structBC->GetInsideCellIndex(i, j, k, is, js, ks, 1);
                    RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);
                    structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                    nxs = lr * xfn(in, jn, kn, iSurface);
                    nys = lr * yfn(in, jn, kn, iSurface);
                    nzs = lr * zfn(in, jn, kn, iSurface);

                    //! calculate the temperature of the inside cell.
                    RDouble staticTemperature = q(is, js, ks, IP) / q(is, js, ks, IR) / coefficientOfStateEquation;

                    RDouble tb = staticTemperature + mfTbCoeff * (q(is, js, ks, IU) * q(is, js, ks, IU) +
                        q(is, js, ks, IV) * q(is, js, ks, IV) + q(is, js, ks, IW) * q(is, js, ks, IW));

                    if (tb < TINY) tb = staticTemperature;
                    RDouble pb = q(is, js, ks, IP) * pow(tb / staticTemperature, refGama / (refGama - 1.0));
                    RDouble rb = pb / tb / coefficientOfStateEquation;

                    RDouble ub = q(is, js, ks, IU) * massFluxRatio;
                    RDouble vb = q(is, js, ks, IV) * massFluxRatio;
                    RDouble wb = q(is, js, ks, IW) * massFluxRatio;

                    //! judge the unusual situation, Vn < 0.
                    RDouble Vn = ub * nxs + vb * nys + wb * nzs;
                    if (Vn < 0.0)
                    {
                        Vn = fluxDensity / rb;
                        ub = Vn * xfn[i];
                        vb = Vn * yfn[i];
                        wb = Vn * zfn[i];
                    }

                    // if the outlet is supersonic flow.
                    RDouble cb = sqrt(refGama * coefficientOfStateEquation * tb);
                    if (Vn >= cb)
                    {
                        rb = q(is, js, ks, IR);
                        ub = q(is, js, ks, IU);
                        vb = q(is, js, ks, IV);
                        wb = q(is, js, ks, IW);
                        pb = q(is, js, ks, IP);
                    }

                    //! allocate value to the ghost cell.
                    for (int layer = 1; layer <= nLayer; ++ layer)
                    {
                        structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);
                        q(it, jt, kt, IR) = rb;
                        q(it, jt, kt, IU) = ub;
                        q(it, jt, kt, IV) = vb;
                        q(it, jt, kt, IW) = wb;
                        q(it, jt, kt, IP) = pb;
                    }
                }
            }
        }
    }

    void NSSolverStruct::NozzleInletFlow(Grid *gridIn, StructBC *structBC)
    {
        using namespace IDX;

        Param_NSSolverStruct *parameters = GetControlParameters();
        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
        Data_Param *bcData = structBC->GetBCParamDataBase();

        RDouble inletVelocity, totalPressure, totalTemperature, direction_inlet[3];
        bcData->GetData("velocity", &inletVelocity, PHDOUBLE, 1);
        bcData->GetData("totalPressure", &totalPressure, PHDOUBLE, 1);
        bcData->GetData("totalTemperature", &totalTemperature, PHDOUBLE, 1);
        bcData->GetData("direction_inlet", &direction_inlet, PHDOUBLE, 3);

        //! Dimensional.
        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
        RDouble refVelocity = refMachNumber * refDimensionalSonicSpeed;
        int isWennScheme = parameters->GetWennSchemeFlag();

        inletVelocity = inletVelocity / refVelocity;
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());

        int nEquation = GetNumberOfEquations();
        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        int iSurface = structBC->GetFaceDirection() + 1;
        int lr = structBC->GetFaceLeftOrRightIndex();
        int is, js, ks, it, jt, kt;
        int in, jn, kn;
        RDouble nxs, nys, nzs;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    structBC->GetInsideCellIndex(i, j, k, is, js, ks, 1);
                    RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);

                    structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                    nxs = lr * xfn(in, jn, kn, iSurface);
                    nys = lr * yfn(in, jn, kn, iSurface);
                    nzs = lr * zfn(in, jn, kn, iSurface);

                    for (int layer = 1; layer <= 1; ++ layer)
                    {
                        structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            q(it, jt, kt, m) = primitiveVarFarfield[m];
                        }
                        q(it, jt, kt, IR) = 0.0;
                        q(it, jt, kt, IU) = inletVelocity * direction_inlet[0];
                        q(it, jt, kt, IV) = inletVelocity * direction_inlet[1];
                        q(it, jt, kt, IW) = inletVelocity * direction_inlet[2];
                        q(it, jt, kt, IP) = 0.0;
                    }

                    int it2, jt2, kt2;
                    structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        q(it2, jt2, kt2, m) = q(it, jt, kt, m);
                    }

                    if (isWennScheme == 1)
                    {
                        int it3, jt3, kt3;
                        structBC->GetGhostCellIndex(i, j, k, it3, jt3, kt3, 3);
                        for (int m = 0; m < nEquation; ++m)
                        {
                            q(it3, jt3, kt3, m) = q(it, jt, kt, m);
                        }
                    }
                }
            }
        }
    }

    void NSSolverStruct::NozzleOutFlow(Grid *gridIn, StructBC *structBC)
    {

    }

    void NSSolverStruct::PressureOutletBC(Grid *gridIn, StructBC *structBC)
    {
        using namespace IDX;

        Param_NSSolverStruct *parameters = GetControlParameters();
        Data_Param *bcData = structBC->GetBCParamDataBase();

        RDouble staticPressure;
        bcData->GetData("staticPressure", &staticPressure, PHDOUBLE, 1);

        //! Dimensional.
        RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
        RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble velocity = refMachNumber * refDimensionalSonicSpeed;
        RDouble dynamicPressure = refDimensionalDensity * velocity * velocity;
        RDouble staticPressureCorf = parameters->GetStaticPressureRelaxCorf();
        staticPressure /= dynamicPressure;

        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int nEquation = GetNumberOfEquations();

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
        RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
        RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());


        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);
        int nSurface = structBC->GetFaceDirection() + 1;
        int leftOrRightIndex = structBC->GetFaceLeftOrRightIndex();

        RDouble normalComponetX, normalComponetY, normalComponetZ;

        int is, js, ks, it, jt, kt;
        int nLayer = GetNumberOfGhostCellLayers();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    structBC->GetInsideCellIndex(i, j, k, is, js, ks, 1);

                    RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);

                    for (int layer = 1; layer <= nLayer; ++ layer)
                    {
                        structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
                        {
                            q(it, jt, kt, iEquation) = q(is, js, ks, iEquation);
                        }
                        //
#ifdef USE_HYFLOW_CODE
                        GAS_SPACE::computeOutPressure(q(it, jt, kt, IP), staticPressure, staticPressureCorf);    //! used in injet , modeified  by dms.
                        if(leftOrRightIndex == -1)
                        {
                            normalComponetX = -faceNormalComponentX(is, js, ks, nSurface);
                            normalComponetY = -faceNormalComponentY(is, js, ks, nSurface);
                            normalComponetZ = -faceNormalComponentZ(is, js, ks, nSurface);
                        }
                        else
                        {
                            normalComponetX = faceNormalComponentX(it, jt, kt, nSurface);
                            normalComponetY = faceNormalComponentY(it, jt, kt, nSurface);
                            normalComponetZ = faceNormalComponentZ(it, jt, kt, nSurface);
                        }

                        GAS_SPACE::computeOutSpeed(q(it, jt, kt, IU), q(it, jt, kt, IV), q(it, jt, kt, IW), normalComponetX, normalComponetY, normalComponetZ);   //! modeified  by dms.

#else
                        q(it, jt, kt, IP) = staticPressure;
#endif
                    }
                }
            }
        }
    }

    void NSSolverStruct::ZeroResiduals(Grid *grid)
    {
        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        residual = 0.0;
    }

    void NSSolverStruct::FillField(Grid *gridIn, FieldProxy *fieldProxy, RDouble valueIn)
    {
        RDouble4D &field = fieldProxy->GetField_STR();
        field = valueIn;
    }

    void NSSolverStruct::FillField(Grid *gridIn, FieldProxy *fieldTargetProxy, FieldProxy *fieldSourceProxy)
    {
        RDouble4D &fieldTarget = fieldTargetProxy->GetField_STR();
        RDouble4D &fieldSource = fieldSourceProxy->GetField_STR();

        fieldTarget = fieldSource;
    }

    FieldProxy *NSSolverStruct::CreateFieldProxy(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Param_NSSolverStruct *parameters = GetControlParameters();
        int isWennScheme = parameters->GetWennSchemeFlag();

        int nEquation = GetNumberOfEquations();

        Range I, J, K;
        GetRange(ni, nj, nk, -2, 1, I, J, K);

        if (isWennScheme == 1)
        {
            GetRange(ni, nj, nk, -3, 2, I, J, K);
        }

        Range M(0, nEquation - 1);

        RDouble4D *field = new RDouble4D(I, J, K, M, fortranArray);

        FieldProxy *field_proxy = new FieldProxy();
        field_proxy->SetField_STR(field, true);

        return field_proxy;
    }

    FieldProxy *NSSolverStruct::GetFieldProxy(Grid *gridIn, const string &fieldName)
    {
        RDouble4D *field = reinterpret_cast <RDouble4D *> (gridIn->GetDataPtr(fieldName));

        FieldProxy *field_proxy = new FieldProxy();

        field_proxy->SetField_STR(field);

        return field_proxy;
    }

    void NSSolverStruct::StoreRhsByResidual(Grid *gridIn, FieldProxy *rightHandSideProxy)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();

        rightHandSide = residual;
    }

    void NSSolverStruct::InitResidual(Grid *gridIn, FieldProxy *rightHandSideProxy)
    {
        StructGrid *grid = StructGridCast(gridIn);

        int nEquation = GetNumberOfEquations();

        Range M(0, nEquation - 1);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        int mst = M.first();
        int med = M.last();

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        residual(i, j, k, m) = -rightHandSide(i, j, k, m);
                    }
                }
            }
        }
    }

    void NSSolverStruct::RecoverResidual(Grid *gridIn, FieldProxy *rightHandSideProxy)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();

        if (grid->GetLevel() != 0)
        {
            residual = rightHandSide;
        }
    }

    FieldProxy *NSSolverStruct::GetResidualProxy(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));

        FieldProxy *residualProxy = new FieldProxy();

        residualProxy->SetField_STR(&residual);

        return residualProxy;
    }

    //! Load flow variables stored in grid to q.
    void NSSolverStruct::LoadQ(Grid *grid, FieldProxy *qProxy)
    {
        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &primitiveVarsOfProxy = qProxy->GetField_STR();
        primitiveVarsOfProxy = primitiveVariables;
    }

    void NSSolverStruct::ComputeViscousCoeff(Grid *gridIn)
    {
        Param_NSSolverStruct *parameters = GetControlParameters();
        int nEnergyRecycle = parameters->GetnEnergyRecycle();
        int nChemical = parameters->GetChemicalFlag();
        if (nChemical == 0)
        {
            ComputeViscousCoefficientWithoutChemical(gridIn);
        }
        else
        {
#ifdef USE_HYFLOW_CODE
            //Compute viscosity.
            //GAS_SPACE::ComputeViscousLaminarCoefficient(StructGridCast(gridIn));
            //Compute transport coefficients.
            if (nEnergyRecycle == 0)
            {
                GAS_SPACE::ComputeTransportCoefficients(StructGridCast(gridIn));
            }
            else
            {
                GAS_SPACE::ComputeTransportCoefficientsR(StructGridCast(gridIn));
            }
#endif
        }
    }

    void NSSolverStruct::ComputeViscousCoefficientWithoutChemical(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        Param_NSSolverStruct *parameters = GetControlParameters();

        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));

        RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();

        //! The visl of ghost cell is computed directly.
        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        using namespace IDX;

        int nGasModel = 0;
        if (GlobalDataBase::IsExist("nGasModel", PHINT, 1))
        {
            nGasModel = GlobalDataBase::GetIntParaFromDB("nGasModel");
        }

        if (nGasModel > 0)
        {
            RDouble primitiveVariables[MAX_SPECIES_NUM], temperature, viscosity;
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        temperature = temperatures(i, j, k, ITT);
                        gas->ComputeViscosityByPrimitive(primitiveVariables, temperature, temperature, nonDimensionalSutherlandTemperature, viscosity);
                        viscousLaminar(i, j, k) = viscosity;
                    }
                }
            }
        }
        else
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        RDouble temperature = temperatures(i, j, k, ITT);
                        viscousLaminar(i, j, k) = temperature * sqrt(temperature) * (1.0 + nonDimensionalSutherlandTemperature) / (temperature + nonDimensionalSutherlandTemperature);
                    }
                }
            }
        }
    }

    void NSSolverStruct::ComputeMinTimeStep(Grid *gridIn, RDouble &minDt, RDouble &maxDt)
    {
        StructGrid *grid = StructGridCast(gridIn);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        minDt = LARGE;
        maxDt = -LARGE;

        RDouble3D &vol = *grid->GetCellVolume();
        RDouble3D &dt = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    minDt = MIN(minDt, static_cast <RDouble> (dt(i, j, k) * vol(i, j, k))); //!when dt*vol is minimum, it doesnot mean that dt is minimum.
                    maxDt = MAX(maxDt, static_cast <RDouble> (dt(i, j, k) * vol(i, j, k))); //!when dt*vol is minimum, it doesnot mean that dt is minimum.
                }
            }
        }
    }

    void NSSolverStruct::ReduceMaxTimeStep(Grid *gridIn, RDouble globalMinDt)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        RDouble3D &vol = *grid->GetCellVolume();
        RDouble3D &dt = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));

        //! ktmax is Used for restricting the time step.
        //! When the value is big, the simulation can be faster but may be divergent. When the value is small, the simulation will be more stable but slower.
        //! The suitable value can improve the stability.
        RDouble ktmax = GlobalDataBase::GetDoubleParaFromDB("ktmax");

        if (ktmax > 0.0)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        dt(i, j, k) = MIN(dt(i, j, k), static_cast <RDouble> (ktmax * globalMinDt / vol(i, j, k)));
                    }
                }
            }
        }
        GhostCell3D(dt, ni, nj, nk);
        //grid->GhostCell3DExceptInterface(dt);
        return;
    }

    void NSSolverStruct::LocalTimeStep(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        RDouble cfl = ComputeCFL();

        RDouble timedt = 0.0;

        RDouble3D &vol = *grid->GetCellVolume();
        RDouble3D &dt = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));
        RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
        RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));

        int nDim = GetDim();
        timedt = 1.0e+30;

        for (int k = kCellStart; k <= kCellEnd; ++ k)     //! Viscous spectrum radius is not considered in computing time step of original code ,modified by liujian 2020-01-09
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    RDouble rad = 0.0;
                    for (int iSurface = 0; iSurface < nDim; ++ iSurface)
                    {
                        rad += invSpectralRadius(i, j, k, iSurface) + visSpectralRadius(i, j, k, iSurface);
                    }
                    dt(i, j, k) = cfl / rad;
                    timedt = min(timedt, dt(i, j, k) * vol(i, j, k));
                }
            }
        }
    }

    void NSSolverStruct::GlobalTimeStep(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble dtau = GlobalDataBase::GetDoubleParaFromDB("dtau");
        RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

        RDouble3D &vol = *grid->GetCellVolume();
        RDouble3D &dt = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    dt(i, j, k) = dtau / reynoldsReferenceLengthDimensional / vol(i, j, k);
                }
            }
        }

        RDouble dtmin = dtau;
        GlobalDataBase::UpdateData("dtmin", &dtmin, PHDOUBLE, 1);
    }

    void NSSolverStruct::LocalGlobalTimeStep(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble cfl = ComputeCFL();

        //! Compute the spectrum radius.
        //SpectrumRadius(grid);

        RDouble3D &vol = *grid->GetCellVolume();
        RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
        RDouble3D &dt = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));

        RDouble rad = 0.0;
        RDouble dtmin = 1.0e30;

        int nDim = GetDim();

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    rad = 0.0;
                    for (int iSurface = 0; iSurface < nDim; ++ iSurface)
                    {
                        rad += invSpectralRadius(i, j, k, iSurface);
                    }
                    RDouble rabc = cfl / rad * vol(i, j, k);

                    dtmin = MIN(dtmin, rabc);
                }
            }
        }

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    dt(i, j, k) = dtmin / vol(i, j, k);
                }
            }
        }
        //grid->GhostCell3DExceptInterface(dt);
        //GhostCell3D(dt, ni, nj, nk);

        GlobalDataBase::UpdateData("dtmin", &dtmin, PHDOUBLE, 1);
    }

    void NSSolverStruct::SpectrumRadius(Grid *grid)
    {
        int ifLowSpeedPrecon = GetControlParameters()->GetIfLowSpeedPrecon();
        if (ifLowSpeedPrecon == 0)
        {
            SpectrumRadiusInviscid(grid);
        }
        else
        {
            PreconSpectrumRadiusInviscid(grid);
        }

        SpectrumRadiusViscous(grid);

#ifdef REMOVE_CODE
        Param_NSSolverStruct *parameters = GetControlParameters();
        int nChemical = parameters->GetChemicalFlag();
        int nChemicalSource = parameters->GetNChemicalSource();
        int nChemicalRadius = parameters->GetNChemicalRadius();

        if (nChemical == 1 && nChemicalSource && nChemicalRadius == 1)
        {
#ifdef USE_HYFLOW_CODE
            GAS_SPACE::AddChemicalSourceInSpectrumRadius(StructGridCast(grid));
#endif
        }
#endif

    }

    void NSSolverStruct::SpectrumRadiusInviscid(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        int RoeEntropyFixMethod = this->GetControlParameters()->GetRoeEntropyFixMethod();

        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());
        RDouble4D &area = *(grid->GetFaceArea());
        RDouble4D &vgn = *(grid->GetFaceNormalVelocity());

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
        RDouble3D &gama = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));
        RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("rtem"));

        invSpectralRadius = 0.0;

        int nDim = GetDim();

        using namespace IDX;

        for (int iSurface = 1; iSurface <= nDim; ++  iSurface)
        {
            int il1, jl1, kl1;
            grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        int il, jl, kl;
                        grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, iSurface);

                        RDouble rm = q(i, j, k, IR);
                        RDouble um = q(i, j, k, IU);
                        RDouble vm = q(i, j, k, IV);
                        RDouble wm = q(i, j, k, IW);
                        RDouble pm = q(i, j, k, IP);

                        RDouble cm = sqrt(gama(i, j, k) * pm / rm);

                        RDouble nx = half * (xfv(i, j, k, iSurface) + xfv(il, jl, kl, iSurface));
                        RDouble ny = half * (yfv(i, j, k, iSurface) + yfv(il, jl, kl, iSurface));
                        RDouble nz = half * (zfv(i, j, k, iSurface) + zfv(il, jl, kl, iSurface));
                        RDouble ub = half * (vgn(i, j, k, iSurface) * area(i, j, k, iSurface) + vgn(il, jl, kl, iSurface) * area(il, jl, kl, iSurface));
                        RDouble vn = nx * um + ny * vm + nz * wm - ub;
                        RDouble ns = sqrt(nx * nx + ny * ny + nz * nz);

                        invSpectralRadius(i, j, k, iSurface - 1) = ABS(vn) + cm * ns;

                        if (RoeEntropyFixMethod == 2)  //! Entropy fix, add by liujian 2019 From JCP ,1995, Lin H-C
                        {
                            RDouble kp = rtem(i, j, k);
                            RDouble Spec_rad_temp = invSpectralRadius(i, j, k, iSurface - 1);

                            RDouble k1 = 0.25;
                            RDouble k2 = 5.0;
                            RDouble eig_lim = Spec_rad_temp * (k1 + k2 * kp);

                            if (Spec_rad_temp < eig_lim)
                            {
                                invSpectralRadius(i, j, k, iSurface - 1) = 0.5 * (Spec_rad_temp * Spec_rad_temp + eig_lim * eig_lim) / (eig_lim + 1.E-20);
                            }
                        }

                    }
                }
            }
        }


        GhostCell3D(invSpectralRadius, ni, nj, nk, nDim);
    }

    void NSSolverStruct::PreconSpectrumRadiusInviscid(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());
        RDouble4D &area = *(grid->GetFaceArea());
        RDouble4D &vgn = *(grid->GetFaceNormalVelocity());

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
        RDouble3D &gama = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));

        RDouble refMachNumber = GetControlParameters()->GetRefMachNumber();

        invSpectralRadius = 0.0;

        int nDim = GetDim();

        using namespace IDX;

        for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
        {
            int il1, jl1, kl1;
            grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        int il, jl, kl;
                        grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, iSurface);

                        RDouble rm = q(i, j, k, IR);
                        RDouble um = q(i, j, k, IU);
                        RDouble vm = q(i, j, k, IV);
                        RDouble wm = q(i, j, k, IW);
                        RDouble pm = q(i, j, k, IP);

                        RDouble cm = sqrt(gama(i, j, k) * pm / rm);

                        RDouble nx = half * (xfn(i, j, k, iSurface) + xfn(il, jl, kl, iSurface));
                        RDouble ny = half * (yfn(i, j, k, iSurface) + yfn(il, jl, kl, iSurface));
                        RDouble nz = half * (zfn(i, j, k, iSurface) + zfn(il, jl, kl, iSurface));
                        RDouble ub = half * (vgn(i, j, k, iSurface) + vgn(il, jl, kl, iSurface));
                        RDouble areaver = half * (area(i, j, k, iSurface) + area(il, jl, kl, iSurface));

                        RDouble vn = nx * um + ny * vm + nz * wm - ub;
                        RDouble ns = sqrt(nx * nx + ny * ny + nz * nz);

                        RDouble v2 = um * um + vm * vm + wm * wm;
                        RDouble c2 = cm * cm;
                        RDouble mach2 = refMachNumber * refMachNumber;

                        RDouble beta = MIN(MAX(v2 / c2, 3.0 * mach2), 1.0);
                        RDouble vnew = half * vn * (1 + beta);
                        RDouble cnew = half * sqrt(((beta - 1) * vn) * ((beta - 1) * vn) + 4 * beta * c2);

                        invSpectralRadius(i, j, k, iSurface - 1) = areaver * (ABS(vnew) + cnew * ns);
                    }
                }
            }
        }

        //GhostCell3D(invSpectralRadius, ni, nj, nk, nDim);
        //grid->GhostCell3DExceptInterface(invSpectralRadius, nDim);
    }

    void NSSolverStruct::SpectrumRadiusViscous(Grid *gridIn)
    {
        //! Compute the spectrum radius of the viscous terms, modified by LiPeng in 2019.
        StructGrid *grid = StructGridCast(gridIn);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());
        RDouble3D &vol = *(grid->GetCellVolume());

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));

        RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));

        RDouble3D &gama = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));

        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));

        visSpectralRadius = 0.0;
        Param_NSSolverStruct *parameters = GetControlParameters();
        int ifViscous = parameters->GetViscousType();
        if (ifViscous == 0)
        {
            return;
        }

        //! Jiri Blazek. Computational fluid dynamics principles and applications (Third Edition)[M]. Elsevier Butterworth-Heinemann, 2015.

        RDouble csrv = 1.0, csrvist = 1.0;
        csrv = GlobalDataBase::GetDoubleParaFromDB("csrv");

        int nDim = GetDim();
        RDouble refReNumber = parameters->GetRefReNumber();
        RDouble prandtlLaminar = parameters->GetPrandtlLaminar();
        RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();
        const RDouble foth = 4.0 / 3.0;

        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int numberOfSpecies = parameters->GetNumberOfSpecies();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nEquation = GetNumberOfEquations();
        int nEnergyRecycle = parameters->GetnEnergyRecycle();

        int mTV = parameters->GetmTV();
        int mTE = parameters->GetmTE();

        RDouble *primitiveVars = new RDouble[nEquation]();
        RDouble *massFractions = NULL;
        RDouble kTurbCoef = 1.0;

        RDouble4D *heatConductivity = NULL, *massDiffusionCoef = NULL;
        RDouble3D *totalCvtr = NULL, *totalCvv = NULL, *totalCve = NULL;
        if (nChemical > 0 && numberOfSpecies > 0)
        {
            massFractions = new RDouble[numberOfSpecies];
            heatConductivity = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("lambda"));
            massDiffusionCoef = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rhoD"));
            totalCvtr = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCvtr"));
            totalCvv = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCvv"));
            totalCve = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCve"));
        }

        using namespace IDX;

#define CFL3D_METHOD


        for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
        {
            int il1, jl1, kl1;
            grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

            for (int k = kCellStart; k <= kCellEnd; ++  k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        int il, jl, kl;
                        grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, iSurface);

                        RDouble nx = half * (xfv(i, j, k, iSurface) + xfv(il, jl, kl, iSurface));
                        RDouble ny = half * (yfv(i, j, k, iSurface) + yfv(il, jl, kl, iSurface));
                        RDouble nz = half * (zfv(i, j, k, iSurface) + zfv(il, jl, kl, iSurface));

                        RDouble ns2 = nx * nx + ny * ny + nz * nz;

                        //visSpectralRadius(i, j, k, iSurface - 1) = ns2;
                        RDouble rm = q(i, j, k, IR);
                        RDouble vvol = vol(i, j, k);

                        RDouble viscLaminar = viscousLaminar(i, j, k);
                        RDouble viscTurbulence = viscousTurbulence(i, j, k);
                        RDouble gam_ma = gama(i, j, k);

                        RDouble coef;
#ifdef CFL3D_METHOD
                        RDouble coef1 = foth * (viscLaminar + csrvist * viscTurbulence);
                        RDouble coef2 = gam_ma * (viscLaminar / prandtlLaminar + viscTurbulence / prandtlTurbulence);
                        coef = MAX(coef1, coef2);
#else
                        coef = (viscLaminar + csrvist * viscTurbulence);
#endif

                        //! Compute chemical diffusion terms, added by LiPeng on Feb. 19, 2019
                        RDouble coef3 = 0.0;
                        if (nChemical == 1)       //! The diffusion coefficient of the chemical reaction.
                        {
                            for (int ivar = 0; ivar < nEquation; ++ ivar)
                            {
                                primitiveVars[ivar] = q(i, j, k, ivar);
                            }

                            for (int ivar = 0; ivar < numberOfSpecies; ++ ivar)
                            {
                                massFractions[ivar] = primitiveVars[nNSEquation + ivar];
                            }
                            //! Obtain the maximum value of (ro*Ds).
                           //coef3 = gas->GetMaximumSpeciesMassDiffusionCoef(primitiveVars, viscLaminar, viscTurbulence);
                            for (int ivar = 0; ivar < gas->GetnLastSpeciesIndex(); ++ ivar)
                            {
                                coef3 = MAX(coef3, (*massDiffusionCoef)(i, j, k, ivar));
                            }

                            //! Multi-temperature model.
                            if (nTemperatureModel > 1)
                            {
                                //! Obtain the temperature.
                                RDouble Tv = t(i, j, k, mTV);
                                RDouble Te = t(i, j, k, mTE);

                                //! Compute heat conductivity./n
                                RDouble Ktr, Kv, Ke;
                                //  gas->ComputeMixtureGasHeatConductivity(primitiveVars, Ttr, Tv, Te, Ktr, Kv, Ke);
                                Ktr = (*heatConductivity)(i, j, k, ITT);
                                Kv = (*heatConductivity)(i, j, k, ITV);
                                RDouble cvtr, cvv, cve;
                                if (nEnergyRecycle == 0)
                                {
                                    cvtr = gas->GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(massFractions);
                                    cvv = gas->GetMixedGasVibrationSpecificHeatAtConstantVolume(primitiveVars, Tv);
                                    cve = gas->GetMixedGasElectronSpecificHeatAtConstantVolume(primitiveVars, Te);
                                }
                                else
                                {
                                    cvtr = (*totalCvtr)(i, j, k);
                                    cvv = (*totalCvv)(i, j, k);
                                    cve = (*totalCve)(i, j, k);
                                }

                                kTurbCoef = 1.0 + viscTurbulence * prandtlLaminar / (viscLaminar * prandtlTurbulence + SMALL);    //modified by dms

                                coef3 = MAX(coef3, Ktr / cvtr * kTurbCoef);
                                RDouble coef4 = 0.0;
                                if (nTemperatureModel == 2)
                                {
                                    Ke = 0.0;
                                    coef4 = (Kv + Ke) / (cvv + cve);
                                }
                                else
                                {
                                    Ke = (*heatConductivity)(i, j, k, ITE);
                                    coef4 = MAX(Kv / cvv, Ke / cve);
                                }
                                coef3 = MAX(coef3, coef4 * kTurbCoef);
                            }
                        }
                        coef = MAX(coef, coef3);
                        //! End computation.

                        coef *= csrv / (refReNumber * rm * vvol);

                        visSpectralRadius(i, j, k, iSurface - 1) = ns2 * coef;

                    }
                }
            }
        }

        delete[]primitiveVars;
        if (nChemical > 0 && numberOfSpecies > 0)
        {
            delete[] massFractions;
        }
    }

    //! Load rhs to residual res stored in grid.
    void NSSolverStruct::LoadResiduals(Grid *gridIn, FieldProxy *rightHandSideProxy)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Range I, J, K;
        GetRange(ni, nj, nk, -2, 1, I, J, K);

        int nEquation = GetNumberOfEquations();
        Range M(0, nEquation - 1);

        RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();
        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));


        residual(I, J, K, M) = -rightHandSide(I, J, K, M);
    }

    RDouble NSSolverStruct::UnsteadyConvergence(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        Param_NSSolverStruct *parameters = GetControlParameters();

        int isUnsteady = parameters->GetIsUnsteady();

        if (!isUnsteady) return zero;

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &qn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
        RDouble4D &qn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));

        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

        RDouble refGama = parameters->GetRefGama();

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        RDouble *prim0 = new RDouble[nEquation];
        RDouble *prim1 = new RDouble[nEquation];
        RDouble *prim2 = new RDouble[nEquation];

        RDouble *qcsv0 = new RDouble[nEquation];
        RDouble *qcsv1 = new RDouble[nEquation];
        RDouble *qcsv2 = new RDouble[nEquation];

        RDouble Ttr[3] = { 0.0 }, Tv[3] = { 0.0 }, Te[3] = { 0.0 };

        using namespace GAS_SPACE;

        RDouble sum1 = zero;
        RDouble sum2 = zero;

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prim0[m] = q(i, j, k, m);
                        prim1[m] = qn1(i, j, k, m);
                        prim2[m] = qn2(i, j, k, m);
                    }

                    if (nTemperatureModel > 1)
                    {
                        gas->GetTemperature(prim0, Ttr[0], Tv[0], Te[0]);
                        gas->GetTemperature(prim1, Ttr[1], Tv[1], Te[1]);
                        gas->GetTemperature(prim2, Ttr[2], Tv[2], Te[2]);
                    }

                    //! Primitive2Conservative function does not need to input gama in fact.
                    //! It will be modified later.
                    gas->Primitive2Conservative(prim0, refGama, Tv[0], Te[0], qcsv0);
                    gas->Primitive2Conservative(prim1, refGama, Tv[1], Te[1], qcsv1);
                    gas->Primitive2Conservative(prim2, refGama, Tv[2], Te[2], qcsv2);


                    for (int m = 0; m < nLaminar; ++ m)
                    {
                        //! Now res has been converted to dq. It means that residual is not right hand side term,but dq.
                        RDouble dq_p = residual(i, j, k, m);        //! qn+1,p+1 - qn+1,p.
                        RDouble dq_n = qcsv0[m] - qcsv1[m];         //! qn+1,p+1 - qn.
                        sum1 += dq_p * dq_p;
                        sum2 += dq_n * dq_n;
                    }
                }
            }
        }

        delete[]prim0;
        delete[]prim1;
        delete[]prim2;
        delete[]qcsv0;
        delete[]qcsv1;
        delete[]qcsv2;

        RDouble cvg = sqrt(ABS(sum1 / (sum2 + SMALL)));
        return cvg;
    }

    void NSSolverStruct::UpdateUnsteadyFlow(Grid *gridIn)
    {
        using namespace IDX;

        Param_NSSolverStruct *parameters = GetControlParameters();

        int isUnsteady = parameters->GetIsUnsteady();
        if (!isUnsteady) return;

        StructGrid *grid = StructGridCast(gridIn);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &qn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
        RDouble4D &qn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));

        RDouble4D &resTmp = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_tmp"));
        RDouble4D &resn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n1"));
        RDouble4D &resn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n2"));

        int nEquation = GetNumberOfEquations();

        qn2 = qn1;
        qn1 = q;

        resn2 = resn1;
        resn1 = resTmp; //! Here the current outnstep is over, the value of the stored resTmp should be assigned to resn1 for the next outnstep. It should be noticed that resTmp only contain the inviscid and viscous flux.

        //voln2 = voln1;
        //voln1 = vol  ;

        //! Statistical variables for DES simulation.
        if (IsNeedStatistics())
        {
            int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
            RDouble statisticalTimePeriod = GlobalDataBase::GetDoubleParaFromDB("statisticalTimePeriod");
            RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

            RDouble c1 = 1.0 / (nStatisticalStep * 1.0);
            if (statisticalTimePeriod > 0.0)
            {
                c1 = MAX(c1, physicalTimeStep / statisticalTimePeriod);
            }
            RDouble c2 = 1.0 - c1;

            RDouble4D &qAverage = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage"));

            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            //qAverage(i, j, k, m) = c1 * qAverage(i, j, k, m) + c2 * q(i, j, k, m);
                            qAverage(i, j, k, m) = c2 * qAverage(i, j, k, m) + c1 * q(i, j, k, m);
                        }
                    }
                }
            }
        }

        //! added by zzp 202108, for computing Reynolds Stress.
        if (IsNeedReynoldsStressStatistics())
        {
            int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
            RDouble statisticalTimePeriod = GlobalDataBase::GetDoubleParaFromDB("statisticalTimePeriod");
            RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

            RDouble c1 = 1.0 / (nStatisticalStep * 1.0);
            if (statisticalTimePeriod > 0.0)
            {
                c1 = MAX(c1, physicalTimeStep / statisticalTimePeriod);
            }
            RDouble c2 = 1.0 - c1;

            RDouble4D &tauAverage = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("tauAverage"));
            RDouble4D &qAverage = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage"));
            RDouble4D &q2Average = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q2Average"));

            int statisticMethod = parameters->GetStatisticMethod();
            if (statisticMethod == 0)
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            RDouble u = q(i, j, k, IU);
                            RDouble v = q(i, j, k, IV);
                            RDouble w = q(i, j, k, IW);

                            q2Average(i, j, k, 0) = c2 * q2Average(i, j, k, 0) + c1 * u * u;
                            q2Average(i, j, k, 1) = c2 * q2Average(i, j, k, 1) + c1 * v * v;
                            q2Average(i, j, k, 2) = c2 * q2Average(i, j, k, 2) + c1 * w * w;
                            q2Average(i, j, k, 3) = c2 * q2Average(i, j, k, 3) + c1 * u * v;
                            q2Average(i, j, k, 4) = c2 * q2Average(i, j, k, 4) + c1 * u * w;
                            q2Average(i, j, k, 5) = c2 * q2Average(i, j, k, 5) + c1 * v * w;
                        }
                    }
                }

                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            RDouble uavg = qAverage(i, j, k, IU);
                            RDouble vavg = qAverage(i, j, k, IV);
                            RDouble wavg = qAverage(i, j, k, IW);

                            tauAverage(i, j, k, 0) = q2Average(i, j, k, 0) - uavg * uavg;
                            tauAverage(i, j, k, 1) = q2Average(i, j, k, 1) - vavg * vavg;
                            tauAverage(i, j, k, 2) = q2Average(i, j, k, 2) - wavg * wavg;
                            tauAverage(i, j, k, 3) = q2Average(i, j, k, 3) - uavg * vavg;
                            tauAverage(i, j, k, 4) = q2Average(i, j, k, 4) - uavg * wavg;
                            tauAverage(i, j, k, 5) = q2Average(i, j, k, 5) - vavg * wavg;
                        }
                    }
                }
            }
            else
            {
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            RDouble uprm = q(i, j, k, IU) - qAverage(i, j, k, IU);
                            RDouble vprm = q(i, j, k, IV) - qAverage(i, j, k, IV);
                            RDouble wprm = q(i, j, k, IW) - qAverage(i, j, k, IW);

                            tauAverage(i, j, k, 0) = c2 * tauAverage(i, j, k, 0) + c1 * uprm * uprm;
                            tauAverage(i, j, k, 1) = c2 * tauAverage(i, j, k, 1) + c1 * vprm * vprm;
                            tauAverage(i, j, k, 2) = c2 * tauAverage(i, j, k, 2) + c1 * wprm * wprm;
                            tauAverage(i, j, k, 3) = c2 * tauAverage(i, j, k, 3) + c1 * uprm * vprm;
                            tauAverage(i, j, k, 4) = c2 * tauAverage(i, j, k, 4) + c1 * uprm * wprm;
                            tauAverage(i, j, k, 5) = c2 * tauAverage(i, j, k, 5) + c1 * vprm * wprm;
                        }
                    }
                }
            }
        }
    }

    void NSSolverStruct::DualTimeSource(Grid *gridIn)
    {
        Param_NSSolverStruct *parameters = GetControlParameters();

        int isUnsteady = parameters->GetIsUnsteady();
        if (!isUnsteady) return;

        StructGrid *grid = StructGridCast(gridIn);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &qn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n1"));
        RDouble4D &qn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_unsteady_n2"));

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &residualn1 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n1"));
        RDouble4D &residualn2 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_n2"));
        RDouble4D &residualTmp = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res_unsteady_tmp"));

        RDouble3D &vol = grid->GetVolume(0); //! The volume at current timestep. If it is dualtime step method, it represents the volume at timestep of n+1.
        RDouble3D &voln1 = grid->GetVolume(1); //! The volume at timestep of n  .
        RDouble3D &voln2 = grid->GetVolume(2); //! the volume at timestep of n-1.

        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

        RDouble refGama = parameters->GetRefGama();

        //! The primitive variables
        RDouble *prim0 = new RDouble[nEquation];
        RDouble *prim1 = new RDouble[nEquation];
        RDouble *prim2 = new RDouble[nEquation];

        //! The conservative variables
        RDouble *qcsv0 = new RDouble[nEquation];
        RDouble *qcsv1 = new RDouble[nEquation];
        RDouble *qcsv2 = new RDouble[nEquation];

        RDouble Ttr[3] = { 0.0 }, Tv[3] = { 0.0 }, Te[3] = { 0.0 };

        using namespace GAS_SPACE;

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        residualTmp = residual;

        //! Computation of dualtime coefficients, including three coefficients for Residual of R(p), R(n) and R(n-1).
        //! and three coefficients for conservative variables of qcsv(p), qcsv(n), and qcsv(n-1).
        RDouble dualTimeCoefficient[7];
        const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);

        RDouble dualTimeResC1 = dualTimeCoefficient[0];
        RDouble dualTimeResC2 = dualTimeCoefficient[1];
        RDouble dualTimeResC3 = dualTimeCoefficient[2];
        RDouble dualTimeQC1 = dualTimeCoefficient[3];
        RDouble dualTimeQC2 = dualTimeCoefficient[4];
        RDouble dualTimeQC3 = dualTimeCoefficient[5];

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prim0[m] = q(i, j, k, m);
                        prim1[m] = qn1(i, j, k, m);
                        prim2[m] = qn2(i, j, k, m);
                    }
                    if (nTemperatureModel > 1)
                    {
                        gas->GetTemperature(prim0, Ttr[0], Tv[0], Te[0]);
                        gas->GetTemperature(prim1, Ttr[1], Tv[1], Te[1]);
                        gas->GetTemperature(prim2, Ttr[2], Tv[2], Te[2]);
                    }
                    //! In fact, there is no need to input the variable of refGama for the function of Primitive2Conservative, it is to be improved.
                    gas->Primitive2Conservative(prim0, refGama, Tv[0], Te[0], qcsv0);
                    gas->Primitive2Conservative(prim1, refGama, Tv[1], Te[1], qcsv1);
                    gas->Primitive2Conservative(prim2, refGama, Tv[2], Te[2], qcsv2);

                    //! Here only for "nLaminar", because the formula is to be improved for "nEquation".
                    for (int m = 0; m < nLaminar; ++ m)
                    {
                        RDouble dualSrcRes = dualTimeResC1 * residual(i, j, k, m) +
                            dualTimeResC2 * residualn1(i, j, k, m) +
                            dualTimeResC3 * residualn2(i, j, k, m);

                        RDouble dualSrcQ = dualTimeQC1 * qcsv0[m] * vol(i, j, k) +
                            dualTimeQC2 * qcsv1[m] * voln1(i, j, k) +
                            dualTimeQC3 * qcsv2[m] * voln2(i, j, k);

                        residual(i, j, k, m) = dualSrcRes + dualSrcQ;
                    }
                }
            }
        }
        delete[]prim0;
        delete[]prim1;
        delete[]prim2;
        delete[]qcsv0;
        delete[]qcsv1;
        delete[]qcsv2;
    }

    void NSSolverStruct::ChemicalSource(Grid *gridIn)
    {
#ifdef USE_HYFLOW_CODE

        StructGrid *grid = StructGridCast(gridIn);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nChemicalRadius = parameters->GetNChemicalRadius();

        if (nChemicalRadius == 1)           //! The source term has been computed in the implicit scheme.
        {
            //GAS_SPACE::ObtainChemicalSource(grid);
            //! Compute chemical source term and its Jacobian matrix.
            int nEnergyRecycle = parameters->GetnEnergyRecycle();
            if (nEnergyRecycle == 0)
            {
                GAS_SPACE::ComputeSourceAndJacobian(grid);
            }
            else
            {
                GAS_SPACE::ComputeSourceAndJacobianD3(grid);
            }

        }
        else                                //! Compute the source term in the explicit scheme.
        {
            GAS_SPACE::ComputeChemicalSource(grid);
        }

#endif 
    }

    //Compute the distance between the center of the grid cell on the first layer and the center of the surface cell.
    void NSSolverStruct::ComputeFirstLayerGridHeight(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
        RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
        RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());
        RDouble2D *firstLayerHeight = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("firstLayerHeight"));
        RDouble4D &faceArea = *(grid->GetFaceArea());
        RDouble3D &gridCellVolume = *(grid->GetCellVolume());

        using namespace IDX;

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        //! Obtain the number of boundary condition regions.
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        int indexOfWall = 0, indexOfCell = 0;

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();

            if (!IsWall(BCType) && BCType != PHENGLEI::ABLATION_SURFACE)
            {
                continue;
            }
            //The following are surface boundary codes.
            int ist, ied, jst, jed, kst, ked;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            indexOfCell = 0;
            int nSurface = structBC->GetFaceDirection() + 1;
            int leftOrRightIndex = structBC->GetFaceLeftOrRightIndex();
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        int iWall, jWall, kWall;
                        structBC->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);

                        //! Obtain the centroid of the surface.
                        RDouble xCenterWall, yCenterWall, zCenterWall;
                        grid->FaceCoor(iWall, jWall, kWall, nSurface, xCenterWall, yCenterWall, zCenterWall);
                        //! Obtain the centroid of the cell on the positive first layer.
                        RDouble xCellCenter, yCellCenter, zCellCenter;
                        grid->CenterCoor(i, j, k, xCellCenter, yCellCenter, zCellCenter);

                        RDouble normalComponetX = leftOrRightIndex * faceNormalComponentX(iWall, jWall, kWall, nSurface);
                        RDouble normalComponetY = leftOrRightIndex * faceNormalComponentY(iWall, jWall, kWall, nSurface);
                        RDouble normalComponetZ = leftOrRightIndex * faceNormalComponentZ(iWall, jWall, kWall, nSurface);

#ifdef USE_ALTERNATIVE_CODE
                        RDouble cellVolume = gridCellVolume(iWall, jWall, kWall); //The volume of first cell.
                        RDouble deltaHeight = half * cellVolume / surfaceArea;
#else
                        //! The projection of the dy called the distance between the center of the cell and that of the surface on normal vector.
                        RDouble deltaHeight = fabs((xCellCenter - xCenterWall) * normalComponetX + (yCellCenter - yCenterWall) * normalComponetY + (zCellCenter - zCenterWall) * normalComponetZ);

                        //RDouble surfaceArea = faceArea(iWall, jWall, kWall, nSurface); // The area of side face.
                        //RDouble cellVolume = gridCellVolume(i, j, k);                  // The volume of first cell.
                        //RDouble deltaHeight = half * cellVolume / surfaceArea;
#endif

                        (*firstLayerHeight)(indexOfWall, indexOfCell) = deltaHeight;
                        ++ indexOfCell; //Next grid cell.
                    }
                }
            }
            ++ indexOfWall; //Next surface regions.
        }

    }

    void NSSolverStruct::ComputeGamaAndTemperature(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &gama = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));

        int nEquation = GetNumberOfEquations();
        Param_NSSolverStruct *parameters = GetControlParameters();
        int nChemical = parameters->GetChemicalFlag();
        int nEquilibriumGas = parameters->GetFlagOfEquilibriumGas();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nSpeciesNumber = parameters->GetNumberOfSpecies();
        int nEnergyRecycle = parameters->GetnEnergyRecycle();

        int mTT = parameters->GetmTT();
        int mTV = parameters->GetmTV();
        int mTE = parameters->GetmTE();

        RDouble refGama = parameters->GetRefGama();

        RDouble4D *speciesEs = NULL, *speciesEnthalpy = NULL, *speciesCvs = NULL, *speciesCps = NULL, *speciesCvvs = NULL, *speciesCves = NULL, *speciesEtrs = NULL, *speciesEvs = NULL, *speciesEes = NULL;
        RDouble3D *totalEnthalpy = NULL, *totalEnergy = NULL, *totalCp = NULL, *totalCv = NULL;
        Thermo_Energy *Thermo_Energy_temparay = NULL;
        RDouble3D *tflag = NULL, *totalCvtr = NULL, *totalCvv = NULL, *totalCve = NULL;
        if (nChemical > 0)
        {
            Thermo_Energy_temparay = gas->GetThermo_Energy_temparay();
            speciesEnthalpy = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("hSpecies"));

            totalEnergy = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalEnergy"));
            totalEnthalpy = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalEnthalpy"));
            totalCv = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCv"));
            totalCp = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCp"));

            tflag = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("tflag"));
            totalCvtr = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCvtr"));
            totalCvv = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCvv"));
            totalCve = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCve"));
            speciesCvs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesCvs"));
            speciesCps = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesCps"));
            speciesCvvs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesCvvs"));
            speciesCves = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesCves"));
            speciesEtrs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEtrs"));
            speciesEvs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEvs"));
            speciesEes = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEes"));
            speciesEs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEs"));
        }

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

        int iInnerCellStart, iInnerCellEnd, jInnerCellStart, jInnerCellEnd, kInnerCellStart, kInnerCellEnd;
        grid->GetCellIterationIndex(iInnerCellStart, iInnerCellEnd, jInnerCellStart, jInnerCellEnd, kInnerCellStart, kInnerCellEnd);

        using namespace GAS_SPACE;
        using namespace IDX;

        if (nChemical == 1)
        {
            RDouble *hSpecies = new RDouble[nSpeciesNumber];
            RDouble *primitiveVars = new RDouble[nEquation];
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            primitiveVars[m] = primitiveVariables(i, j, k, m);
                        }
                        if (nEnergyRecycle == 0)
                        {
                            gas->GetSpecificHeatRatioAndTemperatute(primitiveVars, gama(i, j, k), temperatures(i, j, k, mTT), temperatures(i, j, k, mTV), temperatures(i, j, k, mTE));
                            gas->GetEverySpeciesEnthalpy(temperatures(i, j, k, mTT), temperatures(i, j, k, mTV), temperatures(i, j, k, mTE), hSpecies);

                            for (int s = 0; s < nSpeciesNumber; ++ s)
                            {
                                (*speciesEnthalpy)(i, j, k, s) = hSpecies[s];
                            }
                        }
                        else
                        {
                            if ((*tflag)(i, j, k) < 0.5)
                            {
                                gas->GetTemperatureR(primitiveVars, temperatures(i, j, k, mTT), temperatures(i, j, k, mTV), temperatures(i, j, k, mTE));
                            }

                            gas->ComputeMixturegasThermalParameter(&primitiveVars[IP + 1], temperatures(i, j, k, mTT), temperatures(i, j, k, mTV), temperatures(i, j, k, mTE));
                            gama(i, j, k) = Thermo_Energy_temparay->gama;
                            (*totalCv)(i, j, k) = Thermo_Energy_temparay->Cv;
                            (*totalCp)(i, j, k) = Thermo_Energy_temparay->Cp;
                            (*totalEnergy)(i, j, k) = Thermo_Energy_temparay->E;
                            (*totalEnthalpy)(i, j, k) = Thermo_Energy_temparay->H;

                            for (int s = 0; s < nSpeciesNumber; ++ s)
                            {
                                (*speciesEs)(i, j, k, s) = Thermo_Energy_temparay->Es[s];
                                (*speciesEnthalpy)(i, j, k, s) = Thermo_Energy_temparay->Hs[s];
                                (*speciesCvs)(i, j, k, s) = Thermo_Energy_temparay->Cvs[s];
                                (*speciesCps)(i, j, k, s) = Thermo_Energy_temparay->Cps[s];
                            }

                            if (nTemperatureModel > 1)
                            {
                                (*totalCvtr)(i, j, k) = Thermo_Energy_temparay->Cvtr;
                                (*totalCvv)(i, j, k) = Thermo_Energy_temparay->Cvv;
                                (*totalCve)(i, j, k) = Thermo_Energy_temparay->Cve;

                                for (int s = 0; s < nSpeciesNumber; ++ s)
                                {
                                    (*speciesCvvs)(i, j, k, s) = Thermo_Energy_temparay->Cvvs[s];
                                    (*speciesCves)(i, j, k, s) = Thermo_Energy_temparay->Cves[s];
                                    (*speciesEtrs)(i, j, k, s) = Thermo_Energy_temparay->Etrs[s];
                                    (*speciesEvs)(i, j, k, s) = Thermo_Energy_temparay->Evs[s];
                                    (*speciesEes)(i, j, k, s) = Thermo_Energy_temparay->Ees[s];
                                }
                            }
                        }
                    }
                }
            }
            delete[] hSpecies;
            delete[] primitiveVars;
        }
        else
        {
            if (nEquilibriumGas == 0)
            {
                RDouble omav = one;
                RDouble coefficientOfStateEquation = gas->GetCoefficientOfStateEquation();

                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            RDouble &density = primitiveVariables(i, j, k, IR);
                            RDouble &pressure = primitiveVariables(i, j, k, IP);
                            gama(i, j, k) = refGama;
                            temperatures(i, j, k, ITT) = pressure / (coefficientOfStateEquation * density * omav);
                        }
                    }
                }
            }
#ifdef USE_HYFLOW_CODE
            else //Equilibrium Gas
            {
                RDouble coefficientOfStateEquation = gas->GetCoefficientOfStateEquation();
                RDouble refDensity = parameters->GetRefDimensionalDensity();
                RDouble refVelocity = parameters->GetRefDimensionalVelocity();
                RDouble refTemperature = parameters->GetRefDimensionalTemperature();
                RDouble refPressure = refDensity * refVelocity * refVelocity;
                RDouble gama1 = refGama, density = 0.0, pressure = 0.0, temperature = 0.0;
                RDouble massFraction[11];
                for (int k = kCellStart; k <= kCellEnd; ++ k)
                {
                    for (int j = jCellStart; j <= jCellEnd; ++ j)
                    {
                        for (int i = iCellStart; i <= iCellEnd; ++ i)
                        {
                            //The perfect gas.
                            density = primitiveVariables(i, j, k, IR);
                            pressure = primitiveVariables(i, j, k, IP);
                            temperature = pressure / (coefficientOfStateEquation * density);
                            temperatures(i, j, k, ITT) = temperature;

                            //To compute the equilibrium gas.
                            pressure = pressure * refPressure;
                            temperature = temperature * refTemperature;
                            //Modify the specific heat ratio.
                            eqGas->GetEquilibriumGasComponents(pressure, temperature, massFraction, gama1);
                            gama(i, j, k) = gama1;
                        }
                    }
                }
            }
#endif
        }

        // added by clz : 2012-8-3
        // begin
        //int nTemperatureModel = parameters->GetNTemperatureModel();
        FillCornerPoint3D(temperatures, ni, nj, nk, nTemperatureModel);
        FillCornerPoint3D(gama, ni, nj, nk);
        // end
        // added by clz : 2012-8-3
    }

    void NSSolverStruct::RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coef)
    {
        StructGrid *grid = StructGridCast(gridIn);

        int nEquation = GetNumberOfEquations();

        RDouble4D &dq = dqProxy->GetField_STR();

        RDouble4D &res = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble3D &dt = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        Range M(0, nEquation - 1);
        int mst = M.first();
        int med = M.last();

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {
                        dq(i, j, k, m) = dt(i, j, k) * coef * res(i, j, k, m);
                    }
                }
            }
        }
    }

    void NSSolverStruct::SolutionFix(Grid *gridIn, RDouble *prim, int i, int j, int k)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        int nEquation = GetNumberOfEquations();

        using namespace GAS_SPACE;
        using namespace IDX;

        int np = 0;
        int nDim = GetDim();

        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            prim[iEquation] = zero;
        }

        for (int ii = -1; ii <= 1; ++ ii)
        {
            for (int jj = -1; jj <= 1; ++ jj)
            {
                for (int kk = -1; kk <= 1; ++ kk)
                {
                    int kk2d = kk;
                    if (nDim == 2)
                    {
                        kk2d = 0;
                    }
                    if ((ABS(ii) + ABS(jj) + ABS(kk2d)) == 1)
                    {
                        np = np + 1;
                        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
                        {
                            RDouble f = q(i + ii, j + jj, k + kk2d, iEquation);
                            if (iEquation == IR || iEquation == IP)
                            {
                                f = ABS(f);
                            }
                            prim[iEquation] += f;
                        }
                    }
                }
            }
        }

        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            prim[iEquation] /= np;
        }
    }

    void NSSolverStruct::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
    {
        //return;
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &q = qProxy->GetField_STR();
        RDouble4D &dq = dqProxy->GetField_STR();

        Param_NSSolverStruct *parameters = GetControlParameters();

        RDouble4D &qnew = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &gamma = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));

        RDouble3D &Xm = *(grid->GetCellCenterX());
        RDouble3D &Ym = *(grid->GetCellCenterY());
        RDouble3D &Zm = *(grid->GetCellCenterZ());

        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
        //int nSpecies = nLaminar + nChemical - nNSEquation;
        int nElectronIndex = parameters->GetIndexOfElectron();

        int nDensityModify = parameters->GetnDensityModify();
        int nDebug = parameters->GetnDebug();
        int nEnergyRecycle = parameters->GetnEnergyRecycle();
        RDouble densityMin = parameters->GetDensityMin();
        RDouble densityMinFactor = parameters->GetdensityMinFactor();

        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        RDouble temperatureWall = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
        bool isRadiationEquilibrium = false;
        if (fabs(temperatureWall) <= EPSILON)
        {
            isRadiationEquilibrium = true;
        }
        int mTT = parameters->GetmTT();
        int mTV = parameters->GetmTV();
        int mTE = parameters->GetmTE();

        RDouble3D *totalEnergy = NULL;
        if (nChemical > 0)
        {
            totalEnergy = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalEnergy"));
        }

        RDouble gama = 0.0;
        //! Three modes for temperature, including the translational-rotational temperature, vibrational temperature and electron temperature.
        RDouble rm = 0.0, pm = 0.0, temperature[3] = { 0.0 };

        RDouble *prim = new RDouble[nEquation]();
        RDouble *qtry = new RDouble[nEquation]();
        RDouble *dqlocal = new RDouble[nEquation]();

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        using namespace GAS_SPACE;
        using namespace IDX;

        RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
        RDouble roo = primitiveVarFarfield[IR];
        RDouble poo = primitiveVarFarfield[IP];
        RDouble density_limit = 1.0e-12 * roo;
        RDouble pressure_limit = 1.0e-12 * poo;
        int     nNegativeCell = 0;

        RDouble densityRelax = 1.0;
        RDouble *residualVariation = new RDouble[nEquation];
        RDouble *localMaxResidual = NULL;
        int varIndex = 0;
        int iCodeBranch = GlobalDataBase::GetIntParaFromDB("iCodeBranch");
        if (iCodeBranch > 0)
        {
            varIndex = GlobalDataBase::GetIntParaFromDB("nKeyVariableIndex");
            localMaxResidual = reinterpret_cast <RDouble *> (grid->GetDataPtr("maxResidualVariation"));
            //Initialization.
            for (int m = 0; m < nEquation; ++m)
            {
                localMaxResidual[m] = 0.0;
            }
        }
        
        RDouble3D &xcc = *grid->GetCellCenterX();
        RDouble3D &ycc = *grid->GetCellCenterY();
        RDouble3D &zcc = *grid->GetCellCenterZ();
        int nI = 0, nJ = 0, nK = 0;

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prim[m] = q(i, j, k, m);
                    }
                    gama = gamma(i, j, k);
                    //! Obtain the temperature.
                    temperature[ITT] = t(i, j, k, mTT);
                    temperature[ITV] = t(i, j, k, mTV);
                    temperature[ITE] = t(i, j, k, mTE);

                    if (iCodeBranch > 0)
                    {
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            residualVariation[m] = prim[m];
                        }
                        for (int m = 0; m < 3; ++m)
                        {
                            residualVariation[m + 1] = temperature[m];
                        }
                    }

                    //! Tv and Te is the temperature of the previous iteration step.
                    if (nEnergyRecycle == 0)
                    {
                        gas->Primitive2Conservative(prim, gama, t(i, j, k, mTV), t(i, j, k, mTE), qtry);
                    }
                    else
                    {
                        gas->Primitive2ConservativeR(prim, gama, temperature[ITV], temperature[ITE], (*totalEnergy)(i, j, k), qtry);
                    }

                    if (nDensityModify == 1)
                    {
                        if (dq(i, j, k, 0) > 0.0)
                        {
                            qtry[0] += dq(i, j, k, 0);
                        }
                        else
                        {
                            densityRelax = qtry[0];
                            qtry[0] = qtry[0] * MAX(densityMinFactor, exp(dq(i, j, k, 0) / (qtry[0] + density_limit)));    //! high freestream density.
                            qtry[0] = MAX(densityMin, qtry[0]);
                            densityRelax = (densityRelax - qtry[0] + SMALL) / (-dq(i, j, k, 0) + SMALL);
                        }

                        if (dq(i, j, k, 0) > 0.0 || nChemical == 0)
                        {
                            for (int m = 1; m < nEquation; ++ m)
                            {
                                qtry[m] += dq(i, j, k, m);
                            }
                        }
                        else
                        {
#ifdef USE_HYFLOW_CODE
                            GAS_SPACE::densitiModify(dqlocal, nNSEquation, nLaminar, nEquation, densityRelax);
                            for (int m = 1; m < nEquation; ++ m)
                            {
                                qtry[m] += dq(i, j, k, m) * dqlocal[m];    //! used in chemical.
                            }
#else
                            for (int m = 1; m < nEquation; ++ m)
                            {
                                qtry[m] += dq(i, j, k, m);    //! used in chemical.
                            }
#endif
                        }
                    }
                    else
                    {
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            qtry[m] += dq(i, j, k, m);
                        }
                        qtry[0] = ABS(qtry[0]);
                    }

                    //! As returned parameter, temperature saves the temperatures of the next time-advancing step.
                    if (nEnergyRecycle == 0)
                    {
                        gas->Conservative2Primitive(qtry, gama, prim, temperature);
                    }
                    else
                    {
                        gas->Conservative2PrimitiveR(qtry, gama, prim, temperature);
                    }

                    if (nChemical == 1)    //! Modify the value of dq in the index of nl.
                    {
                        dq(i, j, k, nLaminar) = prim[IR] * prim[nLaminar] - q(i, j, k, IR) * q(i, j, k, nLaminar);
                        if (nElectronIndex >= 0)
                        {
                            dq(i, j, k, nLaminar - 1) = prim[IR] * prim[nLaminar - 1] - q(i, j, k, IR) * q(i, j, k, nLaminar - 1);
                        }
                        if (nTemperatureModel == 2)
                        {
                            dq(i, j, k, nLaminar + 1) = prim[IR] * prim[nLaminar + 1] - q(i, j, k, IR) * q(i, j, k, nLaminar + 1);
                        }
                        else if (nTemperatureModel == 3)
                        {
                            dq(i, j, k, nLaminar + 1) = prim[IR] * prim[nLaminar + 1] - q(i, j, k, IR) * q(i, j, k, nLaminar + 1);
                            dq(i, j, k, nLaminar + 2) = prim[IR] * prim[nLaminar + 2] - q(i, j, k, IR) * q(i, j, k, nLaminar + 2);
                        }
                    }
                    prim[IP] = ABS(prim[IP]);

                    rm = prim[IR];
                    pm = prim[IP];

                    if (rm < density_limit || pm < pressure_limit)
                    {
                        nNegativeCell += 1;
                        SolutionFix(grid, prim, i, j, k);

                        if (nDebug != 0)
                        {
                            int iZone = grid->GetZoneID();
                            int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
                            cout << "  i=   " << i << "  j=    " << j << "   k=    " << k << endl;
                            cout << "  X=   " << Xm(i, j, k) << "  Y=    " << Ym(i, j, k) << "   Z=    " << Zm(i, j, k) << endl;
                            cout << "  density=   " << rm << "  pressure=    " << pm << endl;
                            cout << "  iZone=   " << iZone << "  outIterStep=    " << outIterStep << endl;
                            abort();
                        }
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        qnew(i, j, k, m) = prim[m];
                    }

                    if (iCodeBranch > 0)
                    {
                        //The variation of the variable density.
                        residualVariation[IR] = ABS((prim[IR] - residualVariation[IR])) / (residualVariation[IR] + SMALL);
                        //The variation of the variable pressure.
                        residualVariation[IP] = ABS((prim[IP] - residualVariation[IP])) / (residualVariation[IP] + SMALL);

                        //Compute the variation of the variable temperature.
                        for (int m = 0; m < 3; ++m)
                        {
                            residualVariation[m + 1] = ABS((temperature[m] - residualVariation[m + 1])) / (residualVariation[m + 1] + SMALL);
                        }
                        //Compute the variation of the species mass fractions.
                        for (int m = nNSEquation; m < nEquation; ++m)
                        {
                            residualVariation[m] = ABS((prim[m] - residualVariation[m])) / (residualVariation[m] + SMALL);
                        }

                        //if (residualVariation[varIndex] > localMaxResidual[varIndex])
                        //{
                        //    nI = i;
                        //    nJ = j;
                        //    nK = k;
                        //    localMaxResidual[varIndex] = residualVariation[varIndex];
                        //}

                        for (int m = 0; m < nEquation; ++m)
                        {
                            localMaxResidual[m] = MAX(residualVariation[m], localMaxResidual[m]);
                        }
                    }

                    //! Update the temperature values.
                    if (nEnergyRecycle != 0 || nSlipBCModel > 0 || isRadiationEquilibrium)
                    {
                        for (int m = 0; m < nTemperatureModel; ++ m)
                        {
                            t(i, j, k, m) = temperature[m];
                        }
                    }
                }
            }
        }

        // begin
        if (nNegativeCell > 0)
        {
            cout << "      Warning: negative pressure or density appears in " << nNegativeCell << " cells ... " << endl;
            cout << "               level = " << grid->GetLevel() << endl;
        }
        // end

        delete[] prim;
        delete[] qtry;
        delete[] dqlocal;
        delete[] residualVariation;
    }

    void NSSolverStruct::ViscousFlux(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nChemical = parameters->GetChemicalFlag();
        int nGradPrimtiveMethod = parameters->GetnGradPrimtiveMethod();

        int viscousType = parameters->GetViscousType();
        if (viscousType == INVISCID) return;

        FieldProxy *fluxProxy = CreateFieldProxy(grid);

        int nDim = GetDim();

        int iLES = GlobalDataBase::GetIntParaFromDB("iLES");

        string sgsmodel = " ";
        GlobalDataBase::GetData("sgsmodel", &sgsmodel, PHSTRING, 1);

        if (iLES == NOLES_SOLVER || sgsmodel.substr(0, 6) == "dsmOld")
        {
            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
#ifdef USE_PRE_CODE
                if (nChemical == 0)
                {
                    ReconGradVolWeightedwithCorrection(gridIn, iSurface);

                    CompVisFluxnew(grid, fluxProxy, iSurface);
                }
                else
                {
                    ReconGrad(gridIn, iSurface);

                    CompVisFlux(grid, fluxProxy, iSurface);
                    //! To correct the viscous flux on wall.
                    CorrectViscousFlux(grid, fluxProxy, iSurface);
                }
#else
#ifdef USE_HYFLOW_CODE
                if (nGradPrimtiveMethod == 1)
                {
                    GAS_SPACE::ReconGradVolWeightedwithCorrectionR(gridIn, iSurface);
                }
                else
                {
                    ReconGradVolWeightedwithCorrection(gridIn, iSurface);
                }
#else
                ReconGradVolWeightedwithCorrection(gridIn, iSurface);
#endif

                CompVisFluxnew(grid, fluxProxy, iSurface);

                CorrectViscousFlux(grid, fluxProxy, iSurface);
#endif

                //! To compute the fluxes of cells by plusing the fluxes of all the surrounded faces.
                LoadFlux(grid, fluxProxy, iSurface);
            }
        }
        else
        {
            //ReconGradVolWeightedwithCorrectionLES(gridIn);

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                ReconGradVolWeightedwithCorrection(gridIn, iSurface);

                CompVisFluxLES(grid, fluxProxy, iSurface);

                CorrectViscousFlux(grid, fluxProxy, iSurface);

                //! To compute the fluxes of cells by plusing the fluxes of all the surrounded faces.
                LoadFlux(grid, fluxProxy, iSurface);
            }
        }

        delete fluxProxy;
    }

    void NSSolverStruct::ComputeGradient(Grid *gridIn)
    {

    }

    void NSSolverStruct::ComputeGradientCellCenter(Grid *gridIn)
    {
        StructGrid *strgrid = StructGridCast(gridIn);
        int ni = strgrid->GetNI();
        int nj = strgrid->GetNJ();
        int nk = strgrid->GetNK();

        RDouble4D &xfv = *(strgrid->GetFaceVectorX());
        RDouble4D &yfv = *(strgrid->GetFaceVectorY());
        RDouble4D &zfv = *(strgrid->GetFaceVectorZ());
        RDouble3D &vol = *(strgrid->GetCellVolume());

        int ndim = GetDim();
        Param_NSSolverStruct *parameters = GetControlParameters();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nChemical = parameters->GetChemicalFlag();
        int nNSNumber = parameters->GetNSEquationNumber();
        int nSpeciesNumber = parameters->GetNumberOfSpecies();
        int nGradPrimtiveMethod = parameters->GetnGradPrimtiveMethod();
        int viscousType = parameters->GetViscousType();
        if (nGradPrimtiveMethod == 1 && viscousType < 2)return;

        RDouble4D &primitiveVaiables = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("q"));
        RDouble4D &temperature = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("t"));
        RDouble4D &gradUVWTCellCenterX = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterX"));
        RDouble4D &gradUVWTCellCenterY = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterY"));
        RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterZ"));

        gradUVWTCellCenterX = 0.0;
        gradUVWTCellCenterY = 0.0;
        gradUVWTCellCenterZ = 0.0;

        RDouble4D *gradSpeciesCellCenterX = 0, *gradSpeciesCellCenterY = 0, *gradSpeciesCellCenterZ = 0;
        if (nChemical > 0 && nSpeciesNumber > 0)
        {
            gradSpeciesCellCenterX = reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dcdx"));
            gradSpeciesCellCenterY = reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dcdy"));
            gradSpeciesCellCenterZ = reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dcdz"));
            *gradSpeciesCellCenterX = 0.0;
            *gradSpeciesCellCenterY = 0.0;
            *gradSpeciesCellCenterZ = 0.0;
        }


        for (int nsurf = 1; nsurf <= ndim; ++ nsurf)
        {
            int il1, jl1, kl1;
            strgrid->GetNsurfIndex(il1, jl1, kl1, nsurf);

            int ist = 1;
            int ied = ni - 1 + il1;
            int jst = 1;
            int jed = nj - 1 + jl1;
            int kst = 1;
            int ked = nk - 1 + kl1;

            if (ndim == TWO_D) ked = 1;

 
            //! The gradient of u,v,w.
            for (int m = 0; m < ndim; ++ m)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            int il, jl, kl;
                            strgrid->GetLeftCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, nsurf);

                            RDouble phis = primitiveVaiables(i, j, k, m + 1) + primitiveVaiables(il, jl, kl, m + 1);

                            RDouble ddx = phis * xfv(i, j, k, nsurf);
                            RDouble ddy = phis * yfv(i, j, k, nsurf);
                            RDouble ddz = phis * zfv(i, j, k, nsurf);

                            gradUVWTCellCenterX(i, j, k, m) -= ddx;
                            gradUVWTCellCenterY(i, j, k, m) -= ddy;
                            gradUVWTCellCenterZ(i, j, k, m) -= ddz;

                            gradUVWTCellCenterX(il, jl, kl, m) += ddx;
                            gradUVWTCellCenterY(il, jl, kl, m) += ddy;
                            gradUVWTCellCenterZ(il, jl, kl, m) += ddz;
                        }
                    }
                }
            }

            //! The gradient of temperature.
            for (int m = 0; m < nTemperatureModel; ++ m)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            int il, jl, kl;
                            strgrid->GetLeftCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, nsurf);

                            RDouble phis = temperature(i, j, k, m) + temperature(il, jl, kl, m);

                            RDouble ddx = phis * xfv(i, j, k, nsurf);
                            RDouble ddy = phis * yfv(i, j, k, nsurf);
                            RDouble ddz = phis * zfv(i, j, k, nsurf);

                            gradUVWTCellCenterX(i, j, k, m + 3) -= ddx;
                            gradUVWTCellCenterY(i, j, k, m + 3) -= ddy;
                            gradUVWTCellCenterZ(i, j, k, m + 3) -= ddz;

                            gradUVWTCellCenterX(il, jl, kl, m + 3) += ddx;
                            gradUVWTCellCenterY(il, jl, kl, m + 3) += ddy;
                            gradUVWTCellCenterZ(il, jl, kl, m + 3) += ddz;
                        }
                    }
                }
            }

            if (nChemical == 0)
                continue;

            //! The gradient of mass fraction of species.
            for (int s = 0; s < nSpeciesNumber; ++ s)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            int il, jl, kl;
                            strgrid->GetLeftCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, nsurf);

                            RDouble phis = primitiveVaiables(i, j, k, nNSNumber + s) + primitiveVaiables(il, jl, kl, nNSNumber + s);

                            RDouble ddx = phis * xfv(i, j, k, nsurf);
                            RDouble ddy = phis * yfv(i, j, k, nsurf);
                            RDouble ddz = phis * zfv(i, j, k, nsurf);

                            (*gradSpeciesCellCenterX)(i, j, k, s) -= ddx;
                            (*gradSpeciesCellCenterY)(i, j, k, s) -= ddy;
                            (*gradSpeciesCellCenterZ)(i, j, k, s) -= ddz;

                            (*gradSpeciesCellCenterX)(il, jl, kl, s) += ddx;
                            (*gradSpeciesCellCenterY)(il, jl, kl, s) += ddy;
                            (*gradSpeciesCellCenterZ)(il, jl, kl, s) += ddz;
                        }
                    }
                }
            }

        }

        int ist = 1;
        int ied = ni - 1;
        int jst = 1;
        int jed = nj - 1;
        int kst = 1;
        int ked = nk - 1;



        if (ndim == TWO_D) ked = 1;

        int nTotalVariable = nTemperatureModel + 3;
        for (int m = 0; m < nTotalVariable; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        RDouble oov = half / vol(i, j, k);
                        gradUVWTCellCenterX(i, j, k, m) *= oov;
                        gradUVWTCellCenterY(i, j, k, m) *= oov;
                        gradUVWTCellCenterZ(i, j, k, m) *= oov;
                    }
                }
            }
        }

        if (nChemical > 0)
        {
            for (int m = 0; m < nSpeciesNumber; ++ m)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            RDouble oov = half / vol(i, j, k);
                            (*gradSpeciesCellCenterX)(i, j, k, m) *= oov;
                            (*gradSpeciesCellCenterY)(i, j, k, m) *= oov;
                            (*gradSpeciesCellCenterZ)(i, j, k, m) *= oov;
                        }
                    }
                }
            }
        }
    }

    void NSSolverStruct::ComputeGradientCellCenterForLES(Grid *gridIn)
    {
        StructGrid *strgrid = StructGridCast(gridIn);
        int ni = strgrid->GetNI();
        int nj = strgrid->GetNJ();
        int nk = strgrid->GetNK();

        RDouble4D &xfv = *(strgrid->GetFaceVectorX());
        RDouble4D &yfv = *(strgrid->GetFaceVectorY());
        RDouble4D &zfv = *(strgrid->GetFaceVectorZ());
        RDouble3D &vol = *(strgrid->GetCellVolume());

        int ndim = GetDim();
        Param_NSSolverStruct *parameters = GetControlParameters();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nChemical = parameters->GetChemicalFlag();
        int nNSNumber = parameters->GetNSEquationNumber();
        int nSpeciesNumber = parameters->GetNumberOfSpecies();

        RDouble4D &primitiveVaiables = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("q"));
        RDouble4D &temperature = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("t"));
        RDouble4D &gradUVWTCellCenterX = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterX"));
        RDouble4D &gradUVWTCellCenterY = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterY"));
        RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterZ"));

        gradUVWTCellCenterX = 0.0;
        gradUVWTCellCenterY = 0.0;
        gradUVWTCellCenterZ = 0.0;

        RDouble4D *gradSpeciesCellCenterX = 0, *gradSpeciesCellCenterY = 0, *gradSpeciesCellCenterZ = 0;
        if (nChemical > 0 && nSpeciesNumber > 0)
        {
            gradSpeciesCellCenterX = reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dcdx"));
            gradSpeciesCellCenterY = reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dcdy"));
            gradSpeciesCellCenterZ = reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dcdz"));
            *gradSpeciesCellCenterX = 0.0;
            *gradSpeciesCellCenterY = 0.0;
            *gradSpeciesCellCenterZ = 0.0;
        }


        for (int nsurf = 1; nsurf <= ndim; ++ nsurf)
        {
            int il1, jl1, kl1;
            strgrid->GetNsurfIndex(il1, jl1, kl1, nsurf);

            int ist = 1;
            int ied = ni - 1 + il1;
            int jst = 1;
            int jed = nj - 1 + jl1;
            int kst = 1;
            int ked = nk - 1 + kl1;

            if (ndim == TWO_D) ked = 1;

            //! The gradient of u,v,w.
            for (int m = 0; m < ndim; ++ m)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            int il, jl, kl;
                            strgrid->GetLeftCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, nsurf);

                            RDouble phis = primitiveVaiables(i, j, k, m + 1) + primitiveVaiables(il, jl, kl, m + 1);

                            RDouble ddx = phis * xfv(i, j, k, nsurf);
                            RDouble ddy = phis * yfv(i, j, k, nsurf);
                            RDouble ddz = phis * zfv(i, j, k, nsurf);

                            gradUVWTCellCenterX(i, j, k, m) -= ddx;
                            gradUVWTCellCenterY(i, j, k, m) -= ddy;
                            gradUVWTCellCenterZ(i, j, k, m) -= ddz;

                            gradUVWTCellCenterX(il, jl, kl, m) += ddx;
                            gradUVWTCellCenterY(il, jl, kl, m) += ddy;
                            gradUVWTCellCenterZ(il, jl, kl, m) += ddz;
                        }
                    }
                }
            }

            //! The gradient of temperature.
            for (int m = 0; m < nTemperatureModel; ++ m)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            int il, jl, kl;
                            strgrid->GetLeftCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, nsurf);

                            RDouble phis = temperature(i, j, k, m) + temperature(il, jl, kl, m);

                            RDouble ddx = phis * xfv(i, j, k, nsurf);
                            RDouble ddy = phis * yfv(i, j, k, nsurf);
                            RDouble ddz = phis * zfv(i, j, k, nsurf);

                            gradUVWTCellCenterX(i, j, k, m + 3) -= ddx;
                            gradUVWTCellCenterY(i, j, k, m + 3) -= ddy;
                            gradUVWTCellCenterZ(i, j, k, m + 3) -= ddz;

                            gradUVWTCellCenterX(il, jl, kl, m + 3) += ddx;
                            gradUVWTCellCenterY(il, jl, kl, m + 3) += ddy;
                            gradUVWTCellCenterZ(il, jl, kl, m + 3) += ddz;
                        }
                    }
                }
            }

            if (nChemical == 0)
                continue;

            //! The gradient of mass fraction of species.
            for (int s = 0; s < nSpeciesNumber; ++ s)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            int il, jl, kl;
                            strgrid->GetLeftCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, nsurf);

                            RDouble phis = primitiveVaiables(i, j, k, nNSNumber + s) + primitiveVaiables(il, jl, kl, nNSNumber + s);

                            RDouble ddx = phis * xfv(i, j, k, nsurf);
                            RDouble ddy = phis * yfv(i, j, k, nsurf);
                            RDouble ddz = phis * zfv(i, j, k, nsurf);

                            (*gradSpeciesCellCenterX)(i, j, k, s) -= ddx;
                            (*gradSpeciesCellCenterY)(i, j, k, s) -= ddy;
                            (*gradSpeciesCellCenterZ)(i, j, k, s) -= ddz;

                            (*gradSpeciesCellCenterX)(il, jl, kl, s) += ddx;
                            (*gradSpeciesCellCenterY)(il, jl, kl, s) += ddy;
                            (*gradSpeciesCellCenterZ)(il, jl, kl, s) += ddz;
                        }
                    }
                }
            }

        }
        int ist = 1;
        int ied = ni - 1;
        int jst = 1;
        int jed = nj - 1;
        int kst = 1;
        int ked = nk - 1;


        if (ndim == TWO_D) ked = 1;

        int nTotalVariable = nTemperatureModel + 3;
        for (int m = 0; m < nTotalVariable; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        RDouble oov = half / vol(i, j, k);
                        gradUVWTCellCenterX(i, j, k, m) *= oov;
                        gradUVWTCellCenterY(i, j, k, m) *= oov;
                        gradUVWTCellCenterZ(i, j, k, m) *= oov;
                    }
                }
            }
        }

        GhostCell3D(gradUVWTCellCenterX, ni, nj, nk, ndim);
        GhostCell3D(gradUVWTCellCenterY, ni, nj, nk, ndim);
        GhostCell3D(gradUVWTCellCenterZ, ni, nj, nk, ndim);

        if (nChemical > 0)
        {
            for (int m = 0; m < nSpeciesNumber; ++ m)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            RDouble oov = half / vol(i, j, k);
                            (*gradSpeciesCellCenterX)(i, j, k, m) *= oov;
                            (*gradSpeciesCellCenterY)(i, j, k, m) *= oov;
                            (*gradSpeciesCellCenterZ)(i, j, k, m) *= oov;
                        }
                    }
                }
            }
        }

    }

    void NSSolverStruct::ComputeGradientCellCenterOfruvwptOnlyForMixGrid(Grid *gridIn)
    {
        //! only first layer need these data, to transfer.
        //! added by myk 20201024
        StructGrid *strgrid = StructGridCast(gridIn);
        int ni = strgrid->GetNI();
        int nj = strgrid->GetNJ();
        int nk = strgrid->GetNK();

        RDouble4D &xfv = *(strgrid->GetFaceVectorX());
        RDouble4D &yfv = *(strgrid->GetFaceVectorY());
        RDouble4D &zfv = *(strgrid->GetFaceVectorZ());
        RDouble3D &vol = *(strgrid->GetCellVolume());

        int ndim = GetDim();

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("t"));
        RDouble4D &dqdx = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdx_cc_ruvwpt"));
        RDouble4D &dqdy = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdy_cc_ruvwpt"));
        RDouble4D &dqdz = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdz_cc_ruvwpt"));

        dqdx = 0.0;
        dqdy = 0.0;
        dqdz = 0.0;

        int ist = 1;
        int ied = ni - 1;
        int jst = 1;
        int jed = nj - 1;
        int kst = 1;
        int ked = nk - 1;
        if (ndim == TWO_D) ked = 1;

        int ibc[2];
        int jbc[2];
        int kbc[2];
        ibc[0] = ist;
        jbc[0] = jst;
        kbc[0] = kst;
        ibc[1] = ied;
        jbc[1] = jed;
        kbc[1] = ked;

        if (ndim == TWO_D)
        {
            int k = 1;

            RDouble q_IL[6], q_IR[6];
            RDouble q_JL[6], q_JR[6];

            RDouble xfv_IL, xfv_IR;
            RDouble yfv_IL, yfv_IR;
            RDouble zfv_IL, zfv_IR;

            RDouble xfv_JL, xfv_JR;
            RDouble yfv_JL, yfv_JR;
            RDouble zfv_JL, zfv_JR;
            for (int bcindex = 0; bcindex <= 1; ++ bcindex)
            {
                int j = jbc[bcindex];
                for (int i = ist; i <= ied; ++ i)
                {
                    xfv_IL = xfv(i, j, k, 1);
                    yfv_IL = yfv(i, j, k, 1);
                    zfv_IL = zfv(i, j, k, 1);

                    xfv_IR = xfv(i + 1, j, k, 1);
                    yfv_IR = yfv(i + 1, j, k, 1);
                    zfv_IR = zfv(i + 1, j, k, 1);

                    xfv_JL = xfv(i, j, k, 2);
                    yfv_JL = yfv(i, j, k, 2);
                    zfv_JL = zfv(i, j, k, 2);

                    xfv_JR = xfv(i, j + 1, k, 2);
                    yfv_JR = yfv(i, j + 1, k, 2);
                    zfv_JR = zfv(i, j + 1, k, 2);

                    for (int m = 0; m <= 4; ++ m)
                    {
                        q_IL[m] = half * (q(i, j, k, m) + q(i - 1, j, k, m));
                        q_IR[m] = half * (q(i, j, k, m) + q(i + 1, j, k, m));

                        q_JL[m] = half * (q(i, j, k, m) + q(i, j - 1, k, m));
                        q_JR[m] = half * (q(i, j, k, m) + q(i, j + 1, k, m));
                    }

                    q_IL[5] = half * (t(i, j, k, 0) + t(i - 1, j, k, 0));
                    q_IR[5] = half * (t(i, j, k, 0) + t(i + 1, j, k, 0));

                    q_JL[5] = half * (t(i, j, k, 0) + t(i, j - 1, k, 0));
                    q_JR[5] = half * (t(i, j, k, 0) + t(i, j + 1, k, 0));

                    for (int m = 0; m <= 5; ++ m)
                    {
                        dqdx(i, j, k, m) = (-q_IL[m] * xfv_IL + q_IR[m] * xfv_IR) + (-q_JL[m] * xfv_JL + q_JR[m] * xfv_JR);
                        dqdy(i, j, k, m) = (-q_IL[m] * yfv_IL + q_IR[m] * yfv_IR) + (-q_JL[m] * yfv_JL + q_JR[m] * yfv_JR);
                        dqdz(i, j, k, m) = (-q_IL[m] * zfv_IL + q_IR[m] * zfv_IR) + (-q_JL[m] * zfv_JL + q_JR[m] * zfv_JR);
                    }
                }
            }

            for (int bcindex = 0; bcindex <= 1; ++ bcindex)
            {
                int i = ibc[bcindex];
                for (int j = jst; j <= jed; ++ j)
                {
                    xfv_IL = xfv(i, j, k, 1);
                    yfv_IL = yfv(i, j, k, 1);
                    zfv_IL = zfv(i, j, k, 1);

                    xfv_IR = xfv(i + 1, j, k, 1);
                    yfv_IR = yfv(i + 1, j, k, 1);
                    zfv_IR = zfv(i + 1, j, k, 1);

                    xfv_JL = xfv(i, j, k, 2);
                    yfv_JL = yfv(i, j, k, 2);
                    zfv_JL = zfv(i, j, k, 2);

                    xfv_JR = xfv(i, j + 1, k, 2);
                    yfv_JR = yfv(i, j + 1, k, 2);
                    zfv_JR = zfv(i, j + 1, k, 2);

                    for (int m = 0; m <= 4; ++ m)
                    {
                        q_IL[m] = half * (q(i, j, k, m) + q(i - 1, j, k, m));
                        q_IR[m] = half * (q(i, j, k, m) + q(i + 1, j, k, m));

                        q_JL[m] = half * (q(i, j, k, m) + q(i, j - 1, k, m));
                        q_JR[m] = half * (q(i, j, k, m) + q(i, j + 1, k, m));
                    }

                    q_IL[5] = half * (t(i, j, k, 0) + t(i - 1, j, k, 0));
                    q_IR[5] = half * (t(i, j, k, 0) + t(i + 1, j, k, 0));

                    q_JL[5] = half * (t(i, j, k, 0) + t(i, j - 1, k, 0));
                    q_JR[5] = half * (t(i, j, k, 0) + t(i, j + 1, k, 0));

                    for (int m = 0; m <= 5; ++ m)
                    {
                        dqdx(i, j, k, m) = (-q_IL[m] * xfv_IL + q_IR[m] * xfv_IR) + (-q_JL[m] * xfv_JL + q_JR[m] * xfv_JR);
                        dqdy(i, j, k, m) = (-q_IL[m] * yfv_IL + q_IR[m] * yfv_IR) + (-q_JL[m] * yfv_JL + q_JR[m] * yfv_JR);
                        dqdz(i, j, k, m) = (-q_IL[m] * zfv_IL + q_IR[m] * zfv_IR) + (-q_JL[m] * zfv_JL + q_JR[m] * zfv_JR);
                    }
                }
            }
        }
        else
        {
            RDouble q_IL[6], q_IR[6];
            RDouble q_JL[6], q_JR[6];
            RDouble q_KL[6], q_KR[6];

            RDouble xfv_IL, xfv_IR;
            RDouble yfv_IL, yfv_IR;
            RDouble zfv_IL, zfv_IR;

            RDouble xfv_JL, xfv_JR;
            RDouble yfv_JL, yfv_JR;
            RDouble zfv_JL, zfv_JR;

            RDouble xfv_KL, xfv_KR;
            RDouble yfv_KL, yfv_KR;
            RDouble zfv_KL, zfv_KR;

            for (int bcindex = 0; bcindex <= 1; ++ bcindex)
            {
                int k = kbc[bcindex];
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        xfv_IL = xfv(i, j, k, 1);
                        yfv_IL = yfv(i, j, k, 1);
                        zfv_IL = zfv(i, j, k, 1);

                        xfv_IR = xfv(i + 1, j, k, 1);
                        yfv_IR = yfv(i + 1, j, k, 1);
                        zfv_IR = zfv(i + 1, j, k, 1);

                        xfv_JL = xfv(i, j, k, 2);
                        yfv_JL = yfv(i, j, k, 2);
                        zfv_JL = zfv(i, j, k, 2);

                        xfv_JR = xfv(i, j + 1, k, 2);
                        yfv_JR = yfv(i, j + 1, k, 2);
                        zfv_JR = zfv(i, j + 1, k, 2);

                        xfv_KL = xfv(i, j, k, 3);
                        yfv_KL = yfv(i, j, k, 3);
                        zfv_KL = zfv(i, j, k, 3);

                        xfv_KR = xfv(i, j, k + 1, 3);
                        yfv_KR = yfv(i, j, k + 1, 3);
                        zfv_KR = zfv(i, j, k + 1, 3);

                        for (int m = 0; m <= 4; ++ m)
                        {
                            q_IL[m] = half * (q(i, j, k, m) + q(i - 1, j, k, m));
                            q_IR[m] = half * (q(i, j, k, m) + q(i + 1, j, k, m));

                            q_JL[m] = half * (q(i, j, k, m) + q(i, j - 1, k, m));
                            q_JR[m] = half * (q(i, j, k, m) + q(i, j + 1, k, m));

                            q_KL[m] = half * (q(i, j, k, m) + q(i, j, k - 1, m));
                            q_KR[m] = half * (q(i, j, k, m) + q(i, j, k + 1, m));
                        }

                        q_IL[5] = half * (t(i, j, k, 0) + t(i - 1, j, k, 0));
                        q_IR[5] = half * (t(i, j, k, 0) + t(i + 1, j, k, 0));

                        q_JL[5] = half * (t(i, j, k, 0) + t(i, j - 1, k, 0));
                        q_JR[5] = half * (t(i, j, k, 0) + t(i, j + 1, k, 0));

                        q_KL[5] = half * (t(i, j, k, 0) + t(i, j, k - 1, 0));
                        q_KR[5] = half * (t(i, j, k, 0) + t(i, j, k + 1, 0));

                        for (int m = 0; m <= 5; ++ m)
                        {
                            dqdx(i, j, k, m) = (-q_IL[m] * xfv_IL + q_IR[m] * xfv_IR) + (-q_JL[m] * xfv_JL + q_JR[m] * xfv_JR) + (-q_KL[m] * xfv_KL + q_KR[m] * xfv_KR);
                            dqdy(i, j, k, m) = (-q_IL[m] * yfv_IL + q_IR[m] * yfv_IR) + (-q_JL[m] * yfv_JL + q_JR[m] * yfv_JR) + (-q_KL[m] * yfv_KL + q_KR[m] * yfv_KR);
                            dqdz(i, j, k, m) = (-q_IL[m] * zfv_IL + q_IR[m] * zfv_IR) + (-q_JL[m] * zfv_JL + q_JR[m] * zfv_JR) + (-q_KL[m] * zfv_KL + q_KR[m] * zfv_KR);
                        }
                    }
                }
            }

            for (int bcindex = 0; bcindex <= 1; ++ bcindex)
            {
                int j = jbc[bcindex];
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        xfv_IL = xfv(i, j, k, 1);
                        yfv_IL = yfv(i, j, k, 1);
                        zfv_IL = zfv(i, j, k, 1);

                        xfv_IR = xfv(i + 1, j, k, 1);
                        yfv_IR = yfv(i + 1, j, k, 1);
                        zfv_IR = zfv(i + 1, j, k, 1);

                        xfv_JL = xfv(i, j, k, 2);
                        yfv_JL = yfv(i, j, k, 2);
                        zfv_JL = zfv(i, j, k, 2);

                        xfv_JR = xfv(i, j + 1, k, 2);
                        yfv_JR = yfv(i, j + 1, k, 2);
                        zfv_JR = zfv(i, j + 1, k, 2);

                        xfv_KL = xfv(i, j, k, 3);
                        yfv_KL = yfv(i, j, k, 3);
                        zfv_KL = zfv(i, j, k, 3);

                        xfv_KR = xfv(i, j, k + 1, 3);
                        yfv_KR = yfv(i, j, k + 1, 3);
                        zfv_KR = zfv(i, j, k + 1, 3);

                        for (int m = 0; m <= 4; ++ m)
                        {
                            q_IL[m] = half * (q(i, j, k, m) + q(i - 1, j, k, m));
                            q_IR[m] = half * (q(i, j, k, m) + q(i + 1, j, k, m));

                            q_JL[m] = half * (q(i, j, k, m) + q(i, j - 1, k, m));
                            q_JR[m] = half * (q(i, j, k, m) + q(i, j + 1, k, m));

                            q_KL[m] = half * (q(i, j, k, m) + q(i, j, k - 1, m));
                            q_KR[m] = half * (q(i, j, k, m) + q(i, j, k + 1, m));
                        }

                        q_IL[5] = half * (t(i, j, k, 0) + t(i - 1, j, k, 0));
                        q_IR[5] = half * (t(i, j, k, 0) + t(i + 1, j, k, 0));

                        q_JL[5] = half * (t(i, j, k, 0) + t(i, j - 1, k, 0));
                        q_JR[5] = half * (t(i, j, k, 0) + t(i, j + 1, k, 0));

                        q_KL[5] = half * (t(i, j, k, 0) + t(i, j, k - 1, 0));
                        q_KR[5] = half * (t(i, j, k, 0) + t(i, j, k + 1, 0));

                        for (int m = 0; m <= 5; ++ m)
                        {
                            dqdx(i, j, k, m) = (-q_IL[m] * xfv_IL + q_IR[m] * xfv_IR) + (-q_JL[m] * xfv_JL + q_JR[m] * xfv_JR) + (-q_KL[m] * xfv_KL + q_KR[m] * xfv_KR);
                            dqdy(i, j, k, m) = (-q_IL[m] * yfv_IL + q_IR[m] * yfv_IR) + (-q_JL[m] * yfv_JL + q_JR[m] * yfv_JR) + (-q_KL[m] * yfv_KL + q_KR[m] * yfv_KR);
                            dqdz(i, j, k, m) = (-q_IL[m] * zfv_IL + q_IR[m] * zfv_IR) + (-q_JL[m] * zfv_JL + q_JR[m] * zfv_JR) + (-q_KL[m] * zfv_KL + q_KR[m] * zfv_KR);
                        }
                    }
                }
            }

            for (int bcindex = 0; bcindex <= 1; ++ bcindex)
            {
                int i = ibc[bcindex];
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        xfv_IL = xfv(i, j, k, 1);
                        yfv_IL = yfv(i, j, k, 1);
                        zfv_IL = zfv(i, j, k, 1);

                        xfv_IR = xfv(i + 1, j, k, 1);
                        yfv_IR = yfv(i + 1, j, k, 1);
                        zfv_IR = zfv(i + 1, j, k, 1);

                        xfv_JL = xfv(i, j, k, 2);
                        yfv_JL = yfv(i, j, k, 2);
                        zfv_JL = zfv(i, j, k, 2);

                        xfv_JR = xfv(i, j + 1, k, 2);
                        yfv_JR = yfv(i, j + 1, k, 2);
                        zfv_JR = zfv(i, j + 1, k, 2);

                        xfv_KL = xfv(i, j, k, 3);
                        yfv_KL = yfv(i, j, k, 3);
                        zfv_KL = zfv(i, j, k, 3);

                        xfv_KR = xfv(i, j, k + 1, 3);
                        yfv_KR = yfv(i, j, k + 1, 3);
                        zfv_KR = zfv(i, j, k + 1, 3);

                        for (int m = 0; m <= 4; ++ m)
                        {
                            q_IL[m] = half * (q(i, j, k, m) + q(i - 1, j, k, m));
                            q_IR[m] = half * (q(i, j, k, m) + q(i + 1, j, k, m));

                            q_JL[m] = half * (q(i, j, k, m) + q(i, j - 1, k, m));
                            q_JR[m] = half * (q(i, j, k, m) + q(i, j + 1, k, m));

                            q_KL[m] = half * (q(i, j, k, m) + q(i, j, k - 1, m));
                            q_KR[m] = half * (q(i, j, k, m) + q(i, j, k + 1, m));
                        }

                        q_IL[5] = half * (t(i, j, k, 0) + t(i - 1, j, k, 0));
                        q_IR[5] = half * (t(i, j, k, 0) + t(i + 1, j, k, 0));

                        q_JL[5] = half * (t(i, j, k, 0) + t(i, j - 1, k, 0));
                        q_JR[5] = half * (t(i, j, k, 0) + t(i, j + 1, k, 0));

                        q_KL[5] = half * (t(i, j, k, 0) + t(i, j, k - 1, 0));
                        q_KR[5] = half * (t(i, j, k, 0) + t(i, j, k + 1, 0));

                        for (int m = 0; m <= 5; ++ m)
                        {
                            dqdx(i, j, k, m) = (-q_IL[m] * xfv_IL + q_IR[m] * xfv_IR) + (-q_JL[m] * xfv_JL + q_JR[m] * xfv_JR) + (-q_KL[m] * xfv_KL + q_KR[m] * xfv_KR);
                            dqdy(i, j, k, m) = (-q_IL[m] * yfv_IL + q_IR[m] * yfv_IR) + (-q_JL[m] * yfv_JL + q_JR[m] * yfv_JR) + (-q_KL[m] * yfv_KL + q_KR[m] * yfv_KR);
                            dqdz(i, j, k, m) = (-q_IL[m] * zfv_IL + q_IR[m] * zfv_IR) + (-q_JL[m] * zfv_JL + q_JR[m] * zfv_JR) + (-q_KL[m] * zfv_KL + q_KR[m] * zfv_KR);
                        }
                    }
                }
            }
        }

        for (int m = 0; m <= 5; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        dqdx(i, j, k, m) /= vol(i, j, k);
                        dqdy(i, j, k, m) /= vol(i, j, k);
                        dqdz(i, j, k, m) /= vol(i, j, k);
                    }
                }
            }
        }
    }


    void NSSolverStruct::ReconGrad(Grid *gridIn, int nsurf)
    {
        Param_NSSolverStruct *parameters = GetControlParameters();
        //! This subroutine is developed to calculate the gradients of a.
        //! Variable phi at a surfaces in the direction nsurf.
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int mst = 0;
        int med = GetNumberOfEquations() - 1;

        int nTemperatureModel = parameters->GetTemperatureModel();
        int mtst = 0;
        int mted = nTemperatureModel - 1;

        RDouble4D &q_4d = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t_4d = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));

        RDouble4D &dqdx_4d = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
        RDouble4D &dqdy_4d = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
        RDouble4D &dqdz_4d = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

        RDouble4D &dtdx_4d = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceX"));
        RDouble4D &dtdy_4d = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceY"));
        RDouble4D &dtdz_4d = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceZ"));

        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());
        RDouble3D &vol = *(grid->GetCellVolume());

        Range I, J, K;
        GetRange(ni, nj, nk, 0, 1, I, J, K);

        int ndim = GetDim();

        int index[] = { 1,2,3,1,2 };

        int ns1 = nsurf;
        int ns2 = index[nsurf];
        int ns3 = index[nsurf + 1];

        int il1 = 0;
        int il2 = 0;
        int il3 = 0;
        int jl1 = 0;
        int jl2 = 0;
        int jl3 = 0;
        int kl1 = 0;
        int kl2 = 0;
        int kl3 = 0;

        if (nsurf == 1)
        {
            il1 = 1;
            jl2 = 1;
            kl3 = 1;
        }
        else if (nsurf == 2)
        {
            jl1 = 1;
            kl2 = 1;
            il3 = 1;
        }
        else if (nsurf == 3)
        {
            kl1 = 1;
            il2 = 1;
            jl3 = 1;
        }

        Range M(mst, med);

        dqdx_4d = 0.0;
        dqdy_4d = 0.0;
        dqdz_4d = 0.0;

        dtdx_4d = 0.0;
        dtdy_4d = 0.0;
        dtdz_4d = 0.0;

        Range IW, JW, KW;
        GetRange(ni, nj, nk, -2, 1, IW, JW, KW);

        RDouble3D worksx(IW, JW, KW, fortranArray);
        RDouble3D worksy(IW, JW, KW, fortranArray);
        RDouble3D worksz(IW, JW, KW, fortranArray);
        RDouble3D worksx1(IW, JW, KW, fortranArray);
        RDouble3D worksy1(IW, JW, KW, fortranArray);
        RDouble3D worksz1(IW, JW, KW, fortranArray);
        RDouble3D workqm(IW, JW, KW, fortranArray);
        RDouble3D worktm(IW, JW, KW, fortranArray);

        worksx(I, J, K) = xfv(I, J, K, ns1) + xfv(I - il1, J - jl1, K - kl1, ns1);
        worksy(I, J, K) = yfv(I, J, K, ns1) + yfv(I - il1, J - jl1, K - kl1, ns1);
        worksz(I, J, K) = zfv(I, J, K, ns1) + zfv(I - il1, J - jl1, K - kl1, ns1);

        for (int m = mst; m <= med; ++ m)
        {
            dqdx_4d(I, J, K, m) = -worksx(I, J, K) * q_4d(I - il1, J - jl1, K - kl1, m);
            dqdy_4d(I, J, K, m) = -worksy(I, J, K) * q_4d(I - il1, J - jl1, K - kl1, m);
            dqdz_4d(I, J, K, m) = -worksz(I, J, K) * q_4d(I - il1, J - jl1, K - kl1, m);
        }

        for (int m = mst; m <= med; ++ m)
        {
            dqdx_4d(I - il1, J - jl1, K - kl1, m) += worksx(I, J, K) * q_4d(I - il1, J - jl1, K - kl1, m);
            dqdy_4d(I - il1, J - jl1, K - kl1, m) += worksy(I, J, K) * q_4d(I - il1, J - jl1, K - kl1, m);
            dqdz_4d(I - il1, J - jl1, K - kl1, m) += worksz(I, J, K) * q_4d(I - il1, J - jl1, K - kl1, m);
        }

        for (int m = mtst; m <= mted; ++ m)
        {
            dtdx_4d(I, J, K, m) = -worksx(I, J, K) * t_4d(I - il1, J - jl1, K - kl1, m);
            dtdy_4d(I, J, K, m) = -worksy(I, J, K) * t_4d(I - il1, J - jl1, K - kl1, m);
            dtdz_4d(I, J, K, m) = -worksz(I, J, K) * t_4d(I - il1, J - jl1, K - kl1, m);
        }

        for (int m = mtst; m <= mted; ++ m)
        {
            dtdx_4d(I - il1, J - jl1, K - kl1, m) += worksx(I, J, K) * t_4d(I - il1, J - jl1, K - kl1, m);
            dtdy_4d(I - il1, J - jl1, K - kl1, m) += worksy(I, J, K) * t_4d(I - il1, J - jl1, K - kl1, m);
            dtdz_4d(I - il1, J - jl1, K - kl1, m) += worksz(I, J, K) * t_4d(I - il1, J - jl1, K - kl1, m);
        }

        if ((nsurf != 2) || (ndim != TWO_D))
        {
            worksx(I, J, K) = xfv(I, J, K, ns2) + xfv(I - il1, J - jl1, K - kl1, ns2);
            worksy(I, J, K) = yfv(I, J, K, ns2) + yfv(I - il1, J - jl1, K - kl1, ns2);
            worksz(I, J, K) = zfv(I, J, K, ns2) + zfv(I - il1, J - jl1, K - kl1, ns2);

            for (int m = mst; m <= med; ++ m)
            {
                workqm(I, J, K) = fourth * (q_4d(I, J, K, m) + q_4d(I - il1, J - jl1, K - kl1, m) + q_4d(I - il2, J - jl2, K - kl2, m) + q_4d(I - il1 - il2, J - jl1 - jl2, K - kl1 - kl2, m));
                dqdx_4d(I, J, K, m) -= worksx(I, J, K) * workqm(I, J, K);
                dqdy_4d(I, J, K, m) -= worksy(I, J, K) * workqm(I, J, K);
                dqdz_4d(I, J, K, m) -= worksz(I, J, K) * workqm(I, J, K);

                dqdx_4d(I - il2, J - jl2, K - kl2, m) += worksx(I, J, K) * workqm(I, J, K);
                dqdy_4d(I - il2, J - jl2, K - kl2, m) += worksy(I, J, K) * workqm(I, J, K);
                dqdz_4d(I - il2, J - jl2, K - kl2, m) += worksz(I, J, K) * workqm(I, J, K);
            }

            for (int m = mtst; m <= mted; ++ m)
            {
                worktm(I, J, K) = fourth * (t_4d(I, J, K, m) + t_4d(I - il1, J - jl1, K - kl1, m) + t_4d(I - il2, J - jl2, K - kl2, m) + t_4d(I - il1 - il2, J - jl1 - jl2, K - kl1 - kl2, m));
                dtdx_4d(I, J, K, m) -= worksx(I, J, K) * worktm(I, J, K);
                dtdy_4d(I, J, K, m) -= worksy(I, J, K) * worktm(I, J, K);
                dtdz_4d(I, J, K, m) -= worksz(I, J, K) * worktm(I, J, K);

                dtdx_4d(I - il2, J - jl2, K - kl2, m) += worksx(I, J, K) * worktm(I, J, K);
                dtdy_4d(I - il2, J - jl2, K - kl2, m) += worksy(I, J, K) * worktm(I, J, K);
                dtdz_4d(I - il2, J - jl2, K - kl2, m) += worksz(I, J, K) * worktm(I, J, K);
            }
        }

        if ((nsurf != 1) || (ndim != TWO_D))
        {
            worksx(I, J, K) = xfv(I, J, K, ns3) + xfv(I - il1, J - jl1, K - kl1, ns3);
            worksy(I, J, K) = yfv(I, J, K, ns3) + yfv(I - il1, J - jl1, K - kl1, ns3);
            worksz(I, J, K) = zfv(I, J, K, ns3) + zfv(I - il1, J - jl1, K - kl1, ns3);

            for (int m = mst; m <= med; ++ m)
            {
                workqm(I, J, K) = fourth * (q_4d(I, J, K, m) + q_4d(I - il1, J - jl1, K - kl1, m) + q_4d(I - il3, J - jl3, K - kl3, m) + q_4d(I - il1 - il3, J - jl1 - jl3, K - kl1 - kl3, m));
                dqdx_4d(I, J, K, m) -= worksx(I, J, K) * workqm(I, J, K);
                dqdy_4d(I, J, K, m) -= worksy(I, J, K) * workqm(I, J, K);
                dqdz_4d(I, J, K, m) -= worksz(I, J, K) * workqm(I, J, K);

                dqdx_4d(I - il3, J - jl3, K - kl3, m) += worksx(I, J, K) * workqm(I, J, K);
                dqdy_4d(I - il3, J - jl3, K - kl3, m) += worksy(I, J, K) * workqm(I, J, K);
                dqdz_4d(I - il3, J - jl3, K - kl3, m) += worksz(I, J, K) * workqm(I, J, K);
            }

            for (int m = mtst; m <= mted; ++ m)
            {
                worktm(I, J, K) = fourth * (t_4d(I, J, K, m) + t_4d(I - il1, J - jl1, K - kl1, m) + t_4d(I - il3, J - jl3, K - kl3, m) + t_4d(I - il1 - il3, J - jl1 - jl3, K - kl1 - kl3, m));
                dtdx_4d(I, J, K, m) -= worksx(I, J, K) * worktm(I, J, K);
                dtdy_4d(I, J, K, m) -= worksy(I, J, K) * worktm(I, J, K);
                dtdz_4d(I, J, K, m) -= worksz(I, J, K) * worktm(I, J, K);

                dtdx_4d(I - il3, J - jl3, K - kl3, m) += worksx(I, J, K) * worktm(I, J, K);
                dtdy_4d(I - il3, J - jl3, K - kl3, m) += worksy(I, J, K) * worktm(I, J, K);
                dtdz_4d(I - il3, J - jl3, K - kl3, m) += worksz(I, J, K) * worktm(I, J, K);
            }
        }

        Range I0(1, ni);
        Range J0(1, nj);
        Range K0(1, nk);

        workqm(I0, J0, K0) = 1.0 / (vol(I0, J0, K0) + vol(I0 - il1, J0 - jl1, K0 - kl1));

        for (int m = mst; m <= med; ++ m)
        {
            dqdx_4d(I0, J0, K0, m) *= workqm(I0, J0, K0);
            dqdy_4d(I0, J0, K0, m) *= workqm(I0, J0, K0);
            dqdz_4d(I0, J0, K0, m) *= workqm(I0, J0, K0);
        }

        for (int m = mtst; m <= mted; ++ m)
        {
            dtdx_4d(I0, J0, K0, m) *= workqm(I0, J0, K0);
            dtdy_4d(I0, J0, K0, m) *= workqm(I0, J0, K0);
            dtdz_4d(I0, J0, K0, m) *= workqm(I0, J0, K0);
        }
    }

    void NSSolverStruct::ReconGradVolWeightedwithCorrection(Grid *gridIn, int nSurface)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &primitiveVaiables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperature = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));

        RDouble4D &gradPrimtiveVarFaceX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
        RDouble4D &gradPrimtiveVarFaceY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
        RDouble4D &gradPrimtiveVarFaceZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));
        gradPrimtiveVarFaceX = 0.0;
        gradPrimtiveVarFaceY = 0.0;
        gradPrimtiveVarFaceZ = 0.0;

        RDouble4D &gradTemperatureFaceX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceX"));
        RDouble4D &gradTemperatureFaceY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceY"));
        RDouble4D &gradTemperatureFaceZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceZ"));
        gradTemperatureFaceX = 0.0;
        gradTemperatureFaceY = 0.0;
        gradTemperatureFaceZ = 0.0;

        RDouble4D &gradUVWTCellCenterX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
        RDouble4D &gradUVWTCellCenterY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
        RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nChemical = parameters->GetChemicalFlag();
        int nNSNumber = parameters->GetNSEquationNumber();
        int nSpeciesNumber = parameters->GetNumberOfSpecies();
        int nTemperatureModel = parameters->GetTemperatureModel();
        RDouble4D *gradSpeciesCellCenterX = 0, *gradSpeciesCellCenterY = 0, *gradSpeciesCellCenterZ = 0;
        if (nChemical > 0 && nSpeciesNumber)
        {
            gradSpeciesCellCenterX = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dcdx"));
            gradSpeciesCellCenterY = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dcdy"));
            gradSpeciesCellCenterZ = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dcdz"));
        }

        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());
        RDouble3D &vol = *(grid->GetCellVolume());

        int il1, jl1, kl1;
        grid->GetNsurfIndex(il1, jl1, kl1, nSurface);

        Range I(1, ni - 1 + il1);
        Range J(1, nj - 1 + jl1);
        Range K(1, nk - 1 + kl1);
        if (nk == 1)
        {
            K.setRange(1, 1);
        }

        int nIndex = 0;
        //! Compute gradient of velocity.
        for (int m = 1; m <= 3; ++ m)
        {
            nIndex = m - 1;
            gradPrimtiveVarFaceX(I, J, K, m) = vol(I, J, K) * gradUVWTCellCenterX(I, J, K, nIndex) + vol(I - il1, J - jl1, K - kl1) * gradUVWTCellCenterX(I - il1, J - jl1, K - kl1, nIndex);
            gradPrimtiveVarFaceY(I, J, K, m) = vol(I, J, K) * gradUVWTCellCenterY(I, J, K, nIndex) + vol(I - il1, J - jl1, K - kl1) * gradUVWTCellCenterY(I - il1, J - jl1, K - kl1, nIndex);
            gradPrimtiveVarFaceZ(I, J, K, m) = vol(I, J, K) * gradUVWTCellCenterZ(I, J, K, nIndex) + vol(I - il1, J - jl1, K - kl1) * gradUVWTCellCenterZ(I - il1, J - jl1, K - kl1, nIndex);

            gradPrimtiveVarFaceX(I, J, K, m) = gradPrimtiveVarFaceX(I, J, K, m) + (xfv(I, J, K, nSurface) * (primitiveVaiables(I, J, K, m) - primitiveVaiables(I - il1, J - jl1, K - kl1, m)) - 0.5 * xfv(I + il1, J + jl1, K + kl1, nSurface) * (primitiveVaiables(I + il1, J + jl1, K + kl1, m) - primitiveVaiables(I, J, K, m)) - 0.5 * xfv(I - il1, J - jl1, K - kl1, nSurface) * (primitiveVaiables(I - il1, J - jl1, K - kl1, m) - primitiveVaiables(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, m)));
            gradPrimtiveVarFaceY(I, J, K, m) = gradPrimtiveVarFaceY(I, J, K, m) + (yfv(I, J, K, nSurface) * (primitiveVaiables(I, J, K, m) - primitiveVaiables(I - il1, J - jl1, K - kl1, m)) - 0.5 * yfv(I + il1, J + jl1, K + kl1, nSurface) * (primitiveVaiables(I + il1, J + jl1, K + kl1, m) - primitiveVaiables(I, J, K, m)) - 0.5 * yfv(I - il1, J - jl1, K - kl1, nSurface) * (primitiveVaiables(I - il1, J - jl1, K - kl1, m) - primitiveVaiables(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, m)));
            gradPrimtiveVarFaceZ(I, J, K, m) = gradPrimtiveVarFaceZ(I, J, K, m) + (zfv(I, J, K, nSurface) * (primitiveVaiables(I, J, K, m) - primitiveVaiables(I - il1, J - jl1, K - kl1, m)) - 0.5 * zfv(I + il1, J + jl1, K + kl1, nSurface) * (primitiveVaiables(I + il1, J + jl1, K + kl1, m) - primitiveVaiables(I, J, K, m)) - 0.5 * zfv(I - il1, J - jl1, K - kl1, nSurface) * (primitiveVaiables(I - il1, J - jl1, K - kl1, m) - primitiveVaiables(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, m)));

            gradPrimtiveVarFaceX(I, J, K, m) /= (vol(I, J, K) + vol(I - il1, J - jl1, K - kl1));
            gradPrimtiveVarFaceY(I, J, K, m) /= (vol(I, J, K) + vol(I - il1, J - jl1, K - kl1));
            gradPrimtiveVarFaceZ(I, J, K, m) /= (vol(I, J, K) + vol(I - il1, J - jl1, K - kl1));
        }

        //! Compute gradient of temperature.
        for (int m = 0; m < nTemperatureModel; ++ m)
        {
            nIndex = m + 3;
            gradTemperatureFaceX(I, J, K, m) = vol(I, J, K) * gradUVWTCellCenterX(I, J, K, nIndex) + vol(I - il1, J - jl1, K - kl1) * gradUVWTCellCenterX(I - il1, J - jl1, K - kl1, nIndex);
            gradTemperatureFaceY(I, J, K, m) = vol(I, J, K) * gradUVWTCellCenterY(I, J, K, nIndex) + vol(I - il1, J - jl1, K - kl1) * gradUVWTCellCenterY(I - il1, J - jl1, K - kl1, nIndex);
            gradTemperatureFaceZ(I, J, K, m) = vol(I, J, K) * gradUVWTCellCenterZ(I, J, K, nIndex) + vol(I - il1, J - jl1, K - kl1) * gradUVWTCellCenterZ(I - il1, J - jl1, K - kl1, nIndex);

            gradTemperatureFaceX(I, J, K, m) = gradTemperatureFaceX(I, J, K, m) + (xfv(I, J, K, nSurface) * (temperature(I, J, K, m) - temperature(I - il1, J - jl1, K - kl1, m)) - 0.5 * xfv(I + il1, J + jl1, K + kl1, nSurface) * (temperature(I + il1, J + jl1, K + kl1, m) - temperature(I, J, K, m)) - 0.5 * xfv(I - il1, J - jl1, K - kl1, nSurface) * (temperature(I - il1, J - jl1, K - kl1, m) - temperature(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, m)));
            gradTemperatureFaceY(I, J, K, m) = gradTemperatureFaceY(I, J, K, m) + (yfv(I, J, K, nSurface) * (temperature(I, J, K, m) - temperature(I - il1, J - jl1, K - kl1, m)) - 0.5 * yfv(I + il1, J + jl1, K + kl1, nSurface) * (temperature(I + il1, J + jl1, K + kl1, m) - temperature(I, J, K, m)) - 0.5 * yfv(I - il1, J - jl1, K - kl1, nSurface) * (temperature(I - il1, J - jl1, K - kl1, m) - temperature(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, m)));
            gradTemperatureFaceZ(I, J, K, m) = gradTemperatureFaceZ(I, J, K, m) + (zfv(I, J, K, nSurface) * (temperature(I, J, K, m) - temperature(I - il1, J - jl1, K - kl1, m)) - 0.5 * zfv(I + il1, J + jl1, K + kl1, nSurface) * (temperature(I + il1, J + jl1, K + kl1, m) - temperature(I, J, K, m)) - 0.5 * zfv(I - il1, J - jl1, K - kl1, nSurface) * (temperature(I - il1, J - jl1, K - kl1, m) - temperature(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, m)));

            gradTemperatureFaceX(I, J, K, m) /= (vol(I, J, K) + vol(I - il1, J - jl1, K - kl1));
            gradTemperatureFaceY(I, J, K, m) /= (vol(I, J, K) + vol(I - il1, J - jl1, K - kl1));
            gradTemperatureFaceZ(I, J, K, m) /= (vol(I, J, K) + vol(I - il1, J - jl1, K - kl1));
        }

        //! Compute gradient of species.
        if (nChemical > 0)
        {
            for (int m = nNSNumber; m < nNSNumber + nSpeciesNumber; ++ m)
            {
                nIndex = m - nNSNumber;
                gradPrimtiveVarFaceX(I, J, K, m) = vol(I, J, K) * (*gradSpeciesCellCenterX)(I, J, K, nIndex) + vol(I - il1, J - jl1, K - kl1) * (*gradSpeciesCellCenterX)(I - il1, J - jl1, K - kl1, nIndex);
                gradPrimtiveVarFaceY(I, J, K, m) = vol(I, J, K) * (*gradSpeciesCellCenterY)(I, J, K, nIndex) + vol(I - il1, J - jl1, K - kl1) * (*gradSpeciesCellCenterY)(I - il1, J - jl1, K - kl1, nIndex);
                gradPrimtiveVarFaceZ(I, J, K, m) = vol(I, J, K) * (*gradSpeciesCellCenterZ)(I, J, K, nIndex) + vol(I - il1, J - jl1, K - kl1) * (*gradSpeciesCellCenterZ)(I - il1, J - jl1, K - kl1, nIndex);

                gradPrimtiveVarFaceX(I, J, K, m) = gradPrimtiveVarFaceX(I, J, K, m) + (xfv(I, J, K, nSurface) * (primitiveVaiables(I, J, K, m) - primitiveVaiables(I - il1, J - jl1, K - kl1, m)) - 0.5 * xfv(I + il1, J + jl1, K + kl1, nSurface) * (primitiveVaiables(I + il1, J + jl1, K + kl1, m) - primitiveVaiables(I, J, K, m)) - 0.5 * xfv(I - il1, J - jl1, K - kl1, nSurface) * (primitiveVaiables(I - il1, J - jl1, K - kl1, m) - primitiveVaiables(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, m)));
                gradPrimtiveVarFaceY(I, J, K, m) = gradPrimtiveVarFaceY(I, J, K, m) + (yfv(I, J, K, nSurface) * (primitiveVaiables(I, J, K, m) - primitiveVaiables(I - il1, J - jl1, K - kl1, m)) - 0.5 * yfv(I + il1, J + jl1, K + kl1, nSurface) * (primitiveVaiables(I + il1, J + jl1, K + kl1, m) - primitiveVaiables(I, J, K, m)) - 0.5 * yfv(I - il1, J - jl1, K - kl1, nSurface) * (primitiveVaiables(I - il1, J - jl1, K - kl1, m) - primitiveVaiables(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, m)));
                gradPrimtiveVarFaceZ(I, J, K, m) = gradPrimtiveVarFaceZ(I, J, K, m) + (zfv(I, J, K, nSurface) * (primitiveVaiables(I, J, K, m) - primitiveVaiables(I - il1, J - jl1, K - kl1, m)) - 0.5 * zfv(I + il1, J + jl1, K + kl1, nSurface) * (primitiveVaiables(I + il1, J + jl1, K + kl1, m) - primitiveVaiables(I, J, K, m)) - 0.5 * zfv(I - il1, J - jl1, K - kl1, nSurface) * (primitiveVaiables(I - il1, J - jl1, K - kl1, m) - primitiveVaiables(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, m)));

                gradPrimtiveVarFaceX(I, J, K, m) /= (vol(I, J, K) + vol(I - il1, J - jl1, K - kl1));
                gradPrimtiveVarFaceY(I, J, K, m) /= (vol(I, J, K) + vol(I - il1, J - jl1, K - kl1));
                gradPrimtiveVarFaceZ(I, J, K, m) /= (vol(I, J, K) + vol(I - il1, J - jl1, K - kl1));
            }
        }

        GetGradientAtFace_CorrectionAtPhysicalBoundary(grid, nSurface);
    }

    void NSSolverStruct::GetGradientAtFace_CorrectionAtPhysicalBoundary(Grid *gridin, int iSurface)
    {
        StructGrid *grid = StructGridCast(gridin);

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());

        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());

        RDouble3D &vol = *(grid->GetCellVolume());

        RDouble4D &gradPrimtiveVarFaceX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
        RDouble4D &gradPrimtiveVarFaceY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
        RDouble4D &gradPrimtiveVarFaceZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

        RDouble4D &gradTemperatureFaceX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceX"));
        RDouble4D &gradTemperatureFaceY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceY"));
        RDouble4D &gradTemperatureFaceZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceZ"));

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));

        RDouble4D &gradUVWTCellCenterX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
        RDouble4D &gradUVWTCellCenterY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
        RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nChemical = parameters->GetChemicalFlag();
        int nNSNumber = parameters->GetNSEquationNumber();
        int nSpeciesNumber = parameters->GetNumberOfSpecies();
        int nTemperatureModel = parameters->GetTemperatureModel();
        RDouble4D *gradSpeciesCellCenterX = 0, *gradSpeciesCellCenterY = 0, *gradSpeciesCellCenterZ = 0;
        int nViscosityFluxSublevelModified = parameters->GetnViscosityFluxSublevelModified();
        if (nChemical > 0 && nSpeciesNumber)
        {
            gradSpeciesCellCenterX = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dcdx"));
            gradSpeciesCellCenterY = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dcdy"));
            gradSpeciesCellCenterZ = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dcdz"));
        }

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();

        using namespace PHENGLEI;

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *bcregion = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = bcregion->GetBCType();
            int nsurf_bc = bcregion->GetFaceDirection() + 1;

            if (BCType < 0 || nsurf_bc != iSurface)
            {
                continue;
            }

            else if (IsWall(BCType))
            {
                int ist, ied, jst, jed, kst, ked;
                bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

                int il1, jl1, kl1;
                grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

                int leftOrRightIndex = bcregion->GetFaceLeftOrRightIndex();  //left -1; right 1

                RDouble nx, ny, nz, vhalf, f_kc_et_ct;

                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            int ibc1, jbc1, kbc1;
                            bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                            nx = xfv(ibc1, jbc1, kbc1, iSurface);
                            ny = yfv(ibc1, jbc1, kbc1, iSurface);
                            nz = zfv(ibc1, jbc1, kbc1, iSurface);

                            //! Modify the gradient of velocity.
                            for (int m = 1; m <= 3; ++ m)
                            {
                                f_kc_et_ct = q(ibc1, jbc1, kbc1, m) - q(ibc1 - il1, jbc1 - jl1, kbc1 - kl1, m);
                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, m) = nx * f_kc_et_ct / vol(i, j, k);
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, m) = ny * f_kc_et_ct / vol(i, j, k);
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, m) = nz * f_kc_et_ct / vol(i, j, k);
                            }

                            //! Modify the gradient of temperature.
                            for (int m = 0; m < nTemperatureModel; ++ m)
                            {
                                f_kc_et_ct = t(ibc1, jbc1, kbc1, m) - t(ibc1 - il1, jbc1 - jl1, kbc1 - kl1, m);
                                gradTemperatureFaceX(ibc1, jbc1, kbc1, m) = nx * f_kc_et_ct / vol(i, j, k);
                                gradTemperatureFaceY(ibc1, jbc1, kbc1, m) = ny * f_kc_et_ct / vol(i, j, k);
                                gradTemperatureFaceZ(ibc1, jbc1, kbc1, m) = nz * f_kc_et_ct / vol(i, j, k);
                            }

                            if (nChemical == 0)
                            {
                                continue;
                            }
                            //! Modify the gradient of species.
                            for (int m = nNSNumber; m < nNSNumber + nSpeciesNumber; ++ m)
                            {
                                f_kc_et_ct = q(ibc1, jbc1, kbc1, m) - q(ibc1 - il1, jbc1 - jl1, kbc1 - kl1, m);
                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, m) = nx * f_kc_et_ct / vol(i, j, k);
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, m) = ny * f_kc_et_ct / vol(i, j, k);
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, m) = nz * f_kc_et_ct / vol(i, j, k);
                            }

                            if (nViscosityFluxSublevelModified == 0)
                            {
                                continue;
                            }

                            ibc1 = ibc1 - il1 * leftOrRightIndex;
                            jbc1 = jbc1 - jl1 * leftOrRightIndex;
                            kbc1 = kbc1 - kl1 * leftOrRightIndex;

                            vhalf = 2.0 / (vol(ibc1, jbc1, kbc1) + vol(ibc1 - il1, jbc1 - jl1, kbc1 - kl1));
                            //    vhalf = 1.0 / vol(i,j,k);

                            nx = xfv(ibc1, jbc1, kbc1, iSurface) * vhalf;
                            ny = yfv(ibc1, jbc1, kbc1, iSurface) * vhalf;
                            nz = zfv(ibc1, jbc1, kbc1, iSurface) * vhalf;


                            //! Modify the gradient of velocity.
                            for (int m = 1; m <= 3; ++ m)
                            {
                                f_kc_et_ct = q(ibc1, jbc1, kbc1, m) - q(ibc1 - il1, jbc1 - jl1, kbc1 - kl1, m);
                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, m) = nx * f_kc_et_ct;
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, m) = ny * f_kc_et_ct;
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, m) = nz * f_kc_et_ct;
                            }

                            //! Modify the gradient of temperature.
                            for (int m = 0; m < nTemperatureModel; ++ m)
                            {
                                f_kc_et_ct = t(ibc1, jbc1, kbc1, m) - t(ibc1 - il1, jbc1 - jl1, kbc1 - kl1, m);
                                gradTemperatureFaceX(ibc1, jbc1, kbc1, m) = nx * f_kc_et_ct;
                                gradTemperatureFaceY(ibc1, jbc1, kbc1, m) = ny * f_kc_et_ct;
                                gradTemperatureFaceZ(ibc1, jbc1, kbc1, m) = nz * f_kc_et_ct;
                            }

                            //! Modify the gradient of species.
                            for (int m = nNSNumber; m < nNSNumber + nSpeciesNumber; ++ m)
                            {
                                f_kc_et_ct = q(ibc1, jbc1, kbc1, m) - q(ibc1 - il1, jbc1 - jl1, kbc1 - kl1, m);
                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, m) = nx * f_kc_et_ct;
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, m) = ny * f_kc_et_ct;
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, m) = nz * f_kc_et_ct;
                            }
                        }
                    }
                }
            }

            else if (BCType == SYMMETRY) // symmetry BC
            {
                RDouble xfnSign;
                RDouble yfnSign;
                RDouble zfnSign;

                int ist, ied, jst, jed, kst, ked;
                bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

                xfnSign = ABS(xfn(ist, jst, kst, iSurface)) - 0.1;
                yfnSign = ABS(yfn(ist, jst, kst, iSurface)) - 0.1;
                zfnSign = ABS(zfn(ist, jst, kst, iSurface)) - 0.1;

                if (xfnSign > 0.0) // for Y-Z symmetry plane
                {
                    for (int k = kst; k <= ked; ++ k)
                    {
                        for (int j = jst; j <= jed; ++ j)
                        {
                            for (int i = ist; i <= ied; ++ i)
                            {
                                int ibc1, jbc1, kbc1;
                                bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, 1) = gradUVWTCellCenterX(i, j, k, 0);
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, 1) = 0.0;
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, 1) = 0.0;

                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, 2) = 0.0;
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, 2) = gradUVWTCellCenterY(i, j, k, 1);
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, 2) = gradUVWTCellCenterZ(i, j, k, 1);

                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, 3) = 0.0;
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, 3) = gradUVWTCellCenterY(i, j, k, 2);
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, 3) = gradUVWTCellCenterZ(i, j, k, 2);

                                //! Modify the gradient of temperature.
                                for (int it = 0; it < nTemperatureModel; ++ it)
                                {
                                    gradTemperatureFaceX(ibc1, jbc1, kbc1, it) = 0.0;
                                    gradTemperatureFaceY(ibc1, jbc1, kbc1, it) = gradUVWTCellCenterY(i, j, k, it + 3);
                                    gradTemperatureFaceZ(ibc1, jbc1, kbc1, it) = gradUVWTCellCenterZ(i, j, k, it + 3);
                                }

                                if (nChemical == 0)
                                {
                                    continue;
                                }
                                //! Modify the gradient of species.
                                for (int s = nNSNumber; s < nNSNumber + nSpeciesNumber; ++ s)
                                {
                                    gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, s) = 0.0;
                                    gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, s) = (*gradSpeciesCellCenterY)(i, j, k, s - nNSNumber);
                                    gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, s) = (*gradSpeciesCellCenterZ)(i, j, k, s - nNSNumber);
                                }

                            }
                        }
                    }
                }

                if (yfnSign > 0.0) // for X-Z symmetry plane
                {
                    for (int k = kst; k <= ked; ++ k)
                    {
                        for (int j = jst; j <= jed; ++ j)
                        {
                            for (int i = ist; i <= ied; ++ i)
                            {
                                int ibc1, jbc1, kbc1;
                                bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, 1) = gradUVWTCellCenterX(i, j, k, 0);
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, 1) = 0.0;
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, 1) = gradUVWTCellCenterZ(i, j, k, 0);

                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, 2) = 0.0;
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, 2) = gradUVWTCellCenterY(i, j, k, 1);
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, 2) = 0.0;

                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, 3) = gradUVWTCellCenterX(i, j, k, 2);
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, 3) = 0.0;
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, 3) = gradUVWTCellCenterZ(i, j, k, 2);

                                //! Modify the gradient of temperature.
                                for (int it = 0; it < nTemperatureModel; ++ it)
                                {
                                    gradTemperatureFaceX(ibc1, jbc1, kbc1, it) = gradUVWTCellCenterX(i, j, k, it + 3);
                                    gradTemperatureFaceY(ibc1, jbc1, kbc1, it) = 0.0;
                                    gradTemperatureFaceZ(ibc1, jbc1, kbc1, it) = gradUVWTCellCenterZ(i, j, k, it + 3);
                                }

                                if (nChemical == 0)
                                {
                                    continue;
                                }
                                //! Modify the gradient of species.
                                for (int s = nNSNumber; s < nNSNumber + nSpeciesNumber; ++ s)
                                {
                                    gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, s) = (*gradSpeciesCellCenterX)(i, j, k, s - nNSNumber);
                                    gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, s) = 0.0;
                                    gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, s) = (*gradSpeciesCellCenterZ)(i, j, k, s - nNSNumber);
                                }

                            }
                        }
                    }
                }

                if (zfnSign > 0.0) // for X-Y symmetry plane
                {
                    for (int k = kst; k <= ked; ++ k)
                    {
                        for (int j = jst; j <= jed; ++ j)
                        {
                            for (int i = ist; i <= ied; ++ i)
                            {
                                int ibc1, jbc1, kbc1;
                                bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, 1) = gradUVWTCellCenterX(i, j, k, 0);
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, 1) = gradUVWTCellCenterY(i, j, k, 0);
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, 1) = 0.0;

                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, 2) = gradUVWTCellCenterX(i, j, k, 1);
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, 2) = gradUVWTCellCenterY(i, j, k, 1);
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, 2) = 0.0;

                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, 3) = 0.0;
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, 3) = 0.0;
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, 3) = gradUVWTCellCenterZ(i, j, k, 2);

                                //! Modify the gradient of temperature.
                                for (int it = 0; it < nTemperatureModel; ++ it)
                                {
                                    gradTemperatureFaceX(ibc1, jbc1, kbc1, it) = gradUVWTCellCenterX(i, j, k, it + 3);
                                    gradTemperatureFaceY(ibc1, jbc1, kbc1, it) = gradUVWTCellCenterY(i, j, k, it + 3);
                                    gradTemperatureFaceZ(ibc1, jbc1, kbc1, it) = 0.0;
                                }

                                if (nChemical == 0)
                                {
                                    continue;
                                }
                                //! Modify the gradient of species.
                                for (int s = nNSNumber; s < nNSNumber + nSpeciesNumber; ++ s)
                                {
                                    gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, s) = (*gradSpeciesCellCenterX)(i, j, k, s - nNSNumber);
                                    gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, s) = (*gradSpeciesCellCenterY)(i, j, k, s - nNSNumber);
                                    gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, s) = 0.0;
                                }

                            }
                        }
                    }
                }
            }
            else //! Other boundary conditions.
            {
                int ist, ied, jst, jed, kst, ked;
                bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            int ibc1, jbc1, kbc1;
                            bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                            for (int m = 1; m <= 3; ++ m)
                            {
                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, m) = gradUVWTCellCenterX(i, j, k, m - 1);
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, m) = gradUVWTCellCenterY(i, j, k, m - 1);
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, m) = gradUVWTCellCenterZ(i, j, k, m - 1);
                            }

                            //! Modify the gradient of temperature.
                            for (int it = 0; it < nTemperatureModel; ++ it)
                            {
                                gradTemperatureFaceX(ibc1, jbc1, kbc1, it) = gradUVWTCellCenterX(i, j, k, it + 3);
                                gradTemperatureFaceY(ibc1, jbc1, kbc1, it) = gradUVWTCellCenterY(i, j, k, it + 3);
                                gradTemperatureFaceZ(ibc1, jbc1, kbc1, it) = gradUVWTCellCenterZ(i, j, k, it + 3);
                            }

                            if (nChemical == 0)
                            {
                                continue;
                            }
                            //! Modify the gradient of species.
                            for (int s = nNSNumber; s < nNSNumber + nSpeciesNumber; ++ s)
                            {
                                gradPrimtiveVarFaceX(ibc1, jbc1, kbc1, s) = (*gradSpeciesCellCenterX)(i, j, k, s - nNSNumber);
                                gradPrimtiveVarFaceY(ibc1, jbc1, kbc1, s) = (*gradSpeciesCellCenterY)(i, j, k, s - nNSNumber);
                                gradPrimtiveVarFaceZ(ibc1, jbc1, kbc1, s) = (*gradSpeciesCellCenterZ)(i, j, k, s - nNSNumber);
                            }

                        }
                    }
                }

            }
        }
    }

    void NSSolverStruct::CompVisFlux(Grid *gridIn, FieldProxy *fluxProxy, int nSurface)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &gradqx = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
        RDouble4D &gradqy = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
        RDouble4D &gradqz = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

        RDouble4D &gradtx = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceX"));
        RDouble4D &gradty = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceY"));
        RDouble4D &gradtz = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceZ"));

        RDouble4D &flux = fluxProxy->GetField_STR();

        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
        int numberOfSpecies = parameters->GetNumberOfSpecies();
        int nElectronIndex = parameters->GetIndexOfElectron();

        RDouble refGama = parameters->GetRefGama();

        RDouble oPrandtlLaminar = parameters->GetoPrandtlLaminar();
        RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble oRefReNumber = parameters->GetoRefReNumber();

        int nrokplus = GlobalDataBase::GetIntParaFromDB("nrokplus");

        RDouble4D *q_turb = 0, *aniss = 0;
        if (nrokplus > 0)
        {
            q_turb = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_turb"));
            aniss = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("aniss"));
        }

        RDouble nx = 0, ny = 0, nz = 0;
        RDouble qx = 0, qy = 0, qz = 0, qvex = 0, qvey = 0, qvez = 0;

        RDouble um = 0, vm = 0, wm = 0, tm = 0, tvm = 0, tem = 0;

        RDouble mul = 0, mut = 0, vis = 0, divv2p3 = 0;
        RDouble txx = 0, tyy = 0, tzz = 0;
        RDouble txy = 0, txz = 0, tyz = 0;

        RDouble dudx = 0, dudy = 0, dudz = 0;
        RDouble dvdx = 0, dvdy = 0, dvdz = 0;
        RDouble dwdx = 0, dwdy = 0, dwdz = 0;
        RDouble dtdx = 0, dtdy = 0, dtdz = 0;

        RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);
        RDouble kcp = 0.0, kv = 0.0, ke = 0.0;      //! The heat conductivity.

        RDouble *qlr = new RDouble[nEquation]();
        RDouble *dfsdx = new RDouble[numberOfSpecies]();
        RDouble *dfsdy = new RDouble[numberOfSpecies]();
        RDouble *dfsdz = new RDouble[numberOfSpecies]();
        RDouble *fvis = new RDouble[nEquation]();
        RDouble *work = new RDouble[nEquation]();

        int iFaceStart = 0, iFaceEnd = 0, jFaceStart = 0, jFaceEnd = 0, kFaceStart = 0, kFaceEnd = 0;
        grid->GetFaceIterationIndex(iFaceStart, iFaceEnd, jFaceStart, jFaceEnd, kFaceStart, kFaceEnd, nSurface);

        int il1, jl1, kl1;
        grid->GetNsurfIndex(il1, jl1, kl1, nSurface);

        using namespace IDX;

        for (int k = kFaceStart; k <= kFaceEnd; ++ k)
        {
            for (int j = jFaceStart; j <= jFaceEnd; ++ j)
            {
                for (int i = iFaceStart; i <= iFaceEnd; ++ i)
                {
                    int il, jl, kl;
                    grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, nSurface);

                    //! Note that if we use iSurface,there must be no omission.
                    nx = xfv(il, jl, kl, nSurface);
                    ny = yfv(il, jl, kl, nSurface);
                    nz = zfv(il, jl, kl, nSurface);

                    um = half * (q(il, jl, kl, IU) + q(i, j, k, IU));
                    vm = half * (q(il, jl, kl, IV) + q(i, j, k, IV));
                    wm = half * (q(il, jl, kl, IW) + q(i, j, k, IW));
                    tm = half * (t(il, jl, kl, ITT) + t(i, j, k, ITT));

                    if (nTemperatureModel == 1)
                    {
                        tvm = tm;  tem = tm;
                    }
                    else if (nTemperatureModel == 2)
                    {
                        tvm = half * (t(il, jl, kl, ITV) + t(i, j, k, ITV));
                        tem = tvm;
                    }
                    else if (nTemperatureModel == 3)
                    {
                        tvm = half * (t(il, jl, kl, ITV) + t(i, j, k, ITV));
                        tem = half * (t(il, jl, kl, ITE) + t(i, j, k, ITE));
                    }

                    mul = half * (viscousLaminar(il, jl, kl) + viscousLaminar(i, j, k));
                    mut = half * (viscousTurbulence(il, jl, kl) + viscousTurbulence(i, j, k));
                    //mut = 0.0;
                    vis = mul + mut;

                    dudx = gradqx(il, jl, kl, IU);
                    dudy = gradqy(il, jl, kl, IU);
                    dudz = gradqz(il, jl, kl, IU);

                    dvdx = gradqx(il, jl, kl, IV);
                    dvdy = gradqy(il, jl, kl, IV);
                    dvdz = gradqz(il, jl, kl, IV);

                    dwdx = gradqx(il, jl, kl, IW);
                    dwdy = gradqy(il, jl, kl, IW);
                    dwdz = gradqz(il, jl, kl, IW);

                    dtdx = gradtx(il, jl, kl, ITT);
                    dtdy = gradty(il, jl, kl, ITT);
                    dtdz = gradtz(il, jl, kl, ITT);

                    qx = 0.0;
                    qy = 0.0;
                    qz = 0.0;

                    if (nChemical == 1)
                    {
                        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
                        {
                            dfsdx[ispecies] = gradqx(il, jl, kl, nNSEquation + ispecies);
                            dfsdy[ispecies] = gradqy(il, jl, kl, nNSEquation + ispecies);
                            dfsdz[ispecies] = gradqz(il, jl, kl, nNSEquation + ispecies);
                        }

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            qlr[m] = half * (q(il, jl, kl, m) + q(i, j, k, m));
                        }

                        gas->ComputeSpeciesMassDiffusionFlux(qlr, dfsdx, dfsdy, dfsdz, tm, tvm, tem, mul, mut, nx, ny, nz, qx, qy, qz, fvis);
                        //! Compute the heat conductivity of multi-component gas.
                        gas->ComputeMixtureGasHeatConductivity(qlr, tm, tvm, tem, kcp, kv, ke);
                    }

                    if (nChemical == 0)    //! Compute the heat conductivity of perfect gas.
                    {
                        kcp = (mul * oPrandtlLaminar + mut * oPrandtlTurbulence) * specificHeatAtConstantPressure;
                    }

                    qx += kcp * dtdx;
                    qy += kcp * dtdy;
                    qz += kcp * dtdz;
                    if (nTemperatureModel == 2)
                    {
                        qvex = (kv + ke) * gradtx(il, jl, kl, ITV);
                        qvey = (kv + ke) * gradty(il, jl, kl, ITV);
                        qvez = (kv + ke) * gradtz(il, jl, kl, ITV);
                        fvis[nLaminar + nChemical] += (nx * qvex + ny * qvey + nz * qvez);
                    }
                    else if (nTemperatureModel == 3)
                    {
                        qvex = kv * gradtx(il, jl, kl, ITV);
                        qvey = kv * gradty(il, jl, kl, ITV);
                        qvez = kv * gradtz(il, jl, kl, ITV);
                        fvis[nLaminar + nChemical] += (nx * qvex + ny * qvey + nz * qvez);
                        qx += qvex;
                        qy += qvey;
                        qz += qvez;

                        qvex = ke * gradtx(il, jl, kl, ITE);
                        qvey = ke * gradty(il, jl, kl, ITE);
                        qvez = ke * gradtz(il, jl, kl, ITE);
                        fvis[nLaminar + nChemical + 1] += (nx * qvex + ny * qvey + nz * qvez);
                    }
                    else
                    {
                        qvex = 0.0;
                        qvey = 0.0;
                        qvez = 0.0;
                    }
                    qx += qvex;
                    qy += qvey;
                    qz += qvez;

                    divv2p3 = two3rd * (dudx + dvdy + dwdz);
                    //! Stress components.
                    txx = vis * (two * dudx - divv2p3);
                    tyy = vis * (two * dvdy - divv2p3);
                    tzz = vis * (two * dwdz - divv2p3);
                    txy = vis * (dudy + dvdx);
                    txz = vis * (dudz + dwdx);
                    tyz = vis * (dvdz + dwdy);

                    fvis[IDX::IR] = 0.0;
                    fvis[IDX::IRU] = nx * txx + ny * txy + nz * txz;
                    fvis[IDX::IRV] = nx * txy + ny * tyy + nz * tyz;
                    fvis[IDX::IRW] = nx * txz + ny * tyz + nz * tzz;
                    fvis[IDX::IRE] = um * fvis[IDX::IRU] + vm * fvis[IDX::IRV] + wm * fvis[IDX::IRW] + (nx * qx + ny * qy + nz * qz);

                    for (int m = 0; m < nLaminar; ++ m)
                    {
                        flux(il, jl, kl, m) = -oRefReNumber * fvis[m];
                    }

                    for (int m = nLaminar + nChemical; m < nEquation; ++ m) // Multi-temperature model.
                    {
                        flux(il, jl, kl, m) = -oRefReNumber * fvis[m];
                    }

                    //! Modify the fluxes of the species whose transport equations need not to be computed.
                    if (nChemical == 1)
                    {
                        flux(il, jl, kl, nLaminar) = 0.0;
                        if (nElectronIndex >= 0)
                        {
                            flux(il, jl, kl, nLaminar - 1) = 0.0;
                        }
                    }
                }
            }
        }

        delete[] fvis;     fvis = NULL;
        delete[] work;     work = NULL;
        delete[] qlr;      qlr = NULL;
        delete[] dfsdx;    dfsdx = NULL;
        delete[] dfsdy;    dfsdy = NULL;
        delete[] dfsdz;    dfsdz = NULL;
    }

    void NSSolverStruct::CorrectViscousFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
    {
        Param_NSSolverStruct *parameters = GetControlParameters();
        int nChemical = parameters->GetChemicalFlag();
        //! Apply the original method, it is no need to modify the flux on wall.
        if (nChemical == 0)
        {
            return;
        }

        RDouble wallTemperature = parameters->GetWallTemperature();
        //! The adiabatic wall applies the primary method.
        if (wallTemperature < 0.0)
        {
            return;
        }

        //! Otherwise, to modify the viscous flux on wall.
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
        RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
        RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());

        RDouble4D &faceArea = *(grid->GetFaceArea());
        RDouble3D &gridCellVolume = *(grid->GetCellVolume());

        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
        RDouble2D *surfaceTemperature = reinterpret_cast <RDouble2D *> (gridIn->GetDataPtr("surfaceTemperature"));
        RDouble4D &faceFluxes = fluxProxy->GetField_STR();
        RDouble oRefReNumber = parameters->GetoRefReNumber();

        //! The number of species.
        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nSpecies = nLaminar + nChemical - nNSEquation;
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

        int viscousType = parameters->GetViscousType();
        RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
        RDouble nondimensionalTemperatureWall = wallTemperature / refDimensionalTemperature;
        RDouble gasConstant = gas->GetUniversalGasConstant();

        //! The type of catalytic wall condition.
        RDouble catalyticCoef = parameters->GetCatalyticCoef();
        RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();
        RDouble localSurfaceTemperature[3] = { nondimensionalTemperatureWall, nondimensionalTemperatureWall, nondimensionalTemperatureWall };

        RDouble *wallVariables = new RDouble[nEquation]();    //! Save the variables on wall.
        RDouble *wallFluxes = new RDouble[nEquation]();
        RDouble *speciesEnthalpy = new RDouble[nSpecies]();
        RDouble *densityDiffusivity = new RDouble[nSpecies]();
        RDouble *vibrationEnergy = new RDouble[nSpecies]();   //! Save the internal energies of vibration.
        RDouble *electronEnergy = new RDouble[nSpecies]();    //! Save the internal energies of electron.
        RDouble3D *surfaceMassFraction = 0;
        RDouble2D *injectionVelocity = 0;
        if (nChemical > 0 && nSpecies > 0)         //! The fully catalytic condition.
        {
            surfaceMassFraction = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
        }
        int isInjection = GlobalDataBase::GetIntParaFromDB("isInjection");
        if (isInjection == 1)
        {
            injectionVelocity = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("injectionVelocity"));
        }

        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        RDouble3D *surfaceSlipVariables = NULL;
        if (nSlipBCModel > 0)
        {
            surfaceSlipVariables = reinterpret_cast <RDouble3D *> (gridIn->GetDataPtr("surfaceSlipVariables"));
        }

        using namespace IDX;
        //RDouble normalComponentX, normalComponentY, normalComponentZ, surfaceArea;
        RDouble transportHeatFlux = 0, speciesHeatFlux = 0;
        RDouble laminarViscosity = 0, turbulentViscosity = 0, viscosity = 0;
        RDouble dudx = 0, dudy = 0, dudz = 0, dvdx = 0, dvdy = 0, dvdz = 0, dwdx = 0, dwdy = 0, dwdz = 0;
        RDouble divv2p3 = 0, txx = 0, tyy = 0, tzz = 0, txy = 0, txz = 0, tyz = 0;

        int indexOfWall = 0, indexOfCell = 0;

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int nBCType = structBC->GetBCType();

            if (!(IsWall(nBCType) && viscousType > INVISCID) && nBCType != PHENGLEI::ABLATION_SURFACE)
            {
                continue;
            }

            int *faceDirectionIndex = structBC->GetFaceDirectionIndex();
            //! 0 stands for i-direction, 1 stands for j-direction, and 2 for k-direction.
            int directionIndex = structBC->GetFaceDirection() + 1;
            //! Judge the face is on the left side or right side in the current zone,-1 stands for left side, and 1 for right side.
            int leftOrRightIndex = structBC->GetFaceLeftOrRightIndex();
            leftOrRightIndex *= -1;

            int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            if (directionIndex != iSurface)
            {
                continue;
            }

            int id = 0, jd = 0, kd = 0;
            GetBCFaceIDX(faceDirectionIndex, id, jd, kd);

            kst = kst + kd;
            ked = ked + kd;
            jst = jst + jd;
            jed = jed + jd;
            ist = ist + id;
            ied = ied + id;

            RDouble velocityXWall = 0.0;
            RDouble velocityYWall = 0.0;
            RDouble velocityZWall = 0.0;

            int nSurface = structBC->GetFaceDirection() + 1;

            indexOfCell = 0;
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        //! The index of cell near the wall.
                        int iCell = i - id;
                        int jCell = j - jd;
                        int kCell = k - kd;

                        //! The index of ghost cell on the first layer.
                        int iL = iCell + faceDirectionIndex[0];
                        int jL = jCell + faceDirectionIndex[1];
                        int kL = kCell + faceDirectionIndex[2];

                        if (isInjection == 1)
                        {
                            int in = 0, jn = 0, kn = 0;
                            structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                            RDouble normalComponetX = leftOrRightIndex * faceNormalComponentX(in, jn, kn, nSurface);
                            RDouble normalComponetY = leftOrRightIndex * faceNormalComponentY(in, jn, kn, nSurface);
                            RDouble normalComponetZ = leftOrRightIndex * faceNormalComponentZ(in, jn, kn, nSurface);

                            RDouble wallVelocity = (*injectionVelocity)(indexOfWall, indexOfCell);

                            velocityXWall = wallVelocity * normalComponetX;
                            velocityYWall = wallVelocity * normalComponetY;
                            velocityZWall = wallVelocity * normalComponetZ;    //Z
                        }

                        if (nSlipBCModel > 0)
                        {
                            velocityXWall = (*surfaceSlipVariables)(indexOfWall, indexOfCell, 3);
                            velocityYWall = (*surfaceSlipVariables)(indexOfWall, indexOfCell, 4);
                            velocityZWall = (*surfaceSlipVariables)(indexOfWall, indexOfCell, 5);
                            for (int m = 0; m < 3; ++ m)
                            {
                                localSurfaceTemperature[m] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m);
                            }
                        }
                        else
                        {
                            //! Obtain the local surface temperatures.Tve=Tv=Te=Tw.
                            localSurfaceTemperature[ITT] = (*surfaceTemperature)(indexOfWall, indexOfCell);
                            localSurfaceTemperature[ITV] = localSurfaceTemperature[ITT];
                            localSurfaceTemperature[ITE] = localSurfaceTemperature[ITT];
                        }

                        //! Compute variables on wall.
                        wallVariables[IU] = velocityXWall;
                        wallVariables[IV] = velocityYWall;
                        wallVariables[IW] = velocityZWall;
                        wallVariables[IP] = primitiveVariables(iCell, jCell, kCell, IP);

                        //! Obtain the mass fractions of species on wall.
                        if (ABS(catalyticCoef) <= EPSILON && nBCType != PHENGLEI::ABLATION_SURFACE) //Fully non-catalytic wall.
                        {
                            for (int m = nNSEquation; m < nNSEquation + nSpecies; ++ m)
                            {
                                wallVariables[m] = primitiveVariables(iCell, jCell, kCell, m);
                            }
                        }
                        else //Fully catalytic wall or finite catalytic wall.
                        {
                            for (int m = nNSEquation; m < nNSEquation + nSpecies; ++ m)
                            {
                                wallVariables[m] = (*surfaceMassFraction)(indexOfWall, indexOfCell, m - nNSEquation);
                            }
                        }

                        if (nSlipBCModel > 0) //Modify the mass fractions on wall.
                        {
                            for (int m = 0; m < nSpecies; ++ m)
                            {
                                wallVariables[nNSEquation + m] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m + 6);
                            }
                        }

                        //! Modify density on wall.
                        RDouble massReciprocal = gas->ComputeMolecularWeightReciprocal(wallVariables);
                        wallVariables[IR] = wallVariables[IP] / (gasConstant * localSurfaceTemperature[ITT] * massReciprocal);

                        //! Obtain the normal vector.
                        RDouble normalComponentX = faceNormalComponentX(i, j, k, iSurface);
                        RDouble normalComponentY = faceNormalComponentY(i, j, k, iSurface);
                        RDouble normalComponentZ = faceNormalComponentZ(i, j, k, iSurface);

                        RDouble surfaceArea = faceArea(i, j, k, iSurface); // The area of side face.
                        RDouble cellVolume = gridCellVolume(iCell, jCell, kCell); // The volume of first cell.
                        RDouble deltaHeight = half * cellVolume / surfaceArea;

                        //! Compute viscosity on wall.
                        gas->ComputeViscosityByPrimitive(wallVariables, localSurfaceTemperature[ITT], localSurfaceTemperature[ITE], nonDimensionalSutherlandTemperature, laminarViscosity); //Tve=Tv=Te=Tw
                        turbulentViscosity = half * (viscousTurbulence(iL, jL, kL) + viscousTurbulence(iCell, jCell, kCell));
                        viscosity = laminarViscosity + turbulentViscosity;

                        //! To obtain the enthalpy of each species.
                        gas->GetEverySpeciesEnthalpy(localSurfaceTemperature[ITT], localSurfaceTemperature[ITV], localSurfaceTemperature[ITE], speciesEnthalpy); //Tve=Tv=Te=Tw

                        //! To obtain the mass diffusion coefficients of each species.
                        gas->GetSpeciesMassDiffusionCoef(wallVariables, laminarViscosity, turbulentViscosity, densityDiffusivity);

                        //! To compute heat conductivity of mixture gas using the Eucken formula and Wassilewa formula.
                        RDouble transRotationHeatCoef, vibrationHeatCoef, electronHeatCoef;
                        gas->ComputeMixtureGasHeatConductivity(wallVariables, localSurfaceTemperature[ITT], localSurfaceTemperature[ITV], localSurfaceTemperature[ITE],
                            transRotationHeatCoef, vibrationHeatCoef, electronHeatCoef);//Tve=Tv=Te=Tw

                        //! Compute the heat flux caused by temperature difference.
                        transportHeatFlux = leftOrRightIndex * transRotationHeatCoef * surfaceArea * ((temperatures(iCell, jCell, kCell, ITT) - localSurfaceTemperature[ITT]) / deltaHeight);
                        if (nTemperatureModel == 2)         //! Two-Temperature Model.
                        {
                            wallFluxes[nLaminar + nChemical] = leftOrRightIndex * (vibrationHeatCoef + electronHeatCoef) * surfaceArea * ((temperatures(iCell, jCell, kCell, ITV) - localSurfaceTemperature[ITV]) / deltaHeight);
                            transportHeatFlux += wallFluxes[nLaminar + nChemical];
                        }
                        else if (nTemperatureModel == 3)    //! Three-Temperature Model.
                        {
                            wallFluxes[nLaminar + nChemical] = leftOrRightIndex * vibrationHeatCoef * surfaceArea * ((temperatures(iCell, jCell, kCell, ITV) - localSurfaceTemperature[ITV]) / deltaHeight);
                            wallFluxes[nLaminar + nChemical + 1] = leftOrRightIndex * electronHeatCoef * surfaceArea * ((temperatures(iCell, jCell, kCell, ITE) - localSurfaceTemperature[ITE]) / deltaHeight);
                            transportHeatFlux += wallFluxes[nLaminar + nChemical] + wallFluxes[nLaminar + nChemical + 1];
                        }

                        //! Compute internal energy of vibration.
                        RDouble peDivideDensity = 0.0;
                        //gas->GetElectronPressure(wallVariables, surfaceTemperature[ITE]) / wallVariables[IR];
                        if (nTemperatureModel > 1)
                        {
                            //! Tve=Tv=Te=Tw.
                            gas->ComputeVibrationEnergy(localSurfaceTemperature[ITV], vibrationEnergy);
                            gas->ComputeElectronEnergy(localSurfaceTemperature[ITE], electronEnergy);
                        }

                        //! Compute the heat flux of species density diffusion.
                        speciesHeatFlux = 0.0;
                        for (int m = 0; m < nSpecies; ++ m)
                        {
                            wallFluxes[nNSEquation + m] = leftOrRightIndex * densityDiffusivity[m] * surfaceArea * ((primitiveVariables(iCell, jCell, kCell, nNSEquation + m) - wallVariables[nNSEquation + m]) / deltaHeight);
                            speciesHeatFlux += speciesEnthalpy[m] * wallFluxes[nNSEquation + m];

                            if (nTemperatureModel == 2)
                            {
                                wallFluxes[nLaminar + nChemical] += (vibrationEnergy[m] + electronEnergy[m] + peDivideDensity) * wallFluxes[nNSEquation + m];    //! modified by dms.
                            }
                            else if (nTemperatureModel == 3)
                            {
                                wallFluxes[nLaminar + nChemical] += vibrationEnergy[m] * wallFluxes[nNSEquation + m];
                                wallFluxes[nLaminar + nChemical + 1] += (electronEnergy[m] + peDivideDensity) * wallFluxes[nNSEquation + m];    //! modified by dms.
                            }
                        }

                        RDouble deltaVelocity = leftOrRightIndex * (primitiveVariables(iCell, jCell, kCell, IU) - wallVariables[IU]);
                        dudx = (deltaVelocity / deltaHeight) * normalComponentX;
                        dudy = (deltaVelocity / deltaHeight) * normalComponentY;
                        dudz = (deltaVelocity / deltaHeight) * normalComponentZ;

                        deltaVelocity = leftOrRightIndex * (primitiveVariables(iCell, jCell, kCell, IV) - wallVariables[IV]);
                        dvdx = (deltaVelocity / deltaHeight) * normalComponentX;
                        dvdy = (deltaVelocity / deltaHeight) * normalComponentY;
                        dvdz = (deltaVelocity / deltaHeight) * normalComponentZ;

                        deltaVelocity = leftOrRightIndex * (primitiveVariables(iCell, jCell, kCell, IW) - wallVariables[IW]);
                        dwdx = (deltaVelocity / deltaHeight) * normalComponentX;
                        dwdy = (deltaVelocity / deltaHeight) * normalComponentY;
                        dwdz = (deltaVelocity / deltaHeight) * normalComponentZ;

                        divv2p3 = 2.0 / 3.0 * (dudx + dvdy + dwdz);
                        //! Compute stress components.
                        txx = viscosity * (2.0 * dudx - divv2p3);
                        tyy = viscosity * (2.0 * dvdy - divv2p3);
                        tzz = viscosity * (2.0 * dwdz - divv2p3);
                        txy = viscosity * (dudy + dvdx);
                        txz = viscosity * (dudz + dwdx);
                        tyz = viscosity * (dvdz + dwdy);

                        wallFluxes[IR] = 0.0;
                        wallFluxes[IRU] = (normalComponentX * txx + normalComponentY * txy + normalComponentZ * txz) * surfaceArea;
                        wallFluxes[IRV] = (normalComponentX * txy + normalComponentY * tyy + normalComponentZ * tyz) * surfaceArea;
                        wallFluxes[IRW] = (normalComponentX * txz + normalComponentY * tyz + normalComponentZ * tzz) * surfaceArea;
                        wallFluxes[IRE] = wallVariables[IU] * wallFluxes[IRU] + wallVariables[IV] * wallFluxes[IRV] + wallVariables[IW] * wallFluxes[IRW] + transportHeatFlux + speciesHeatFlux;

                        for (int m = 0; m < nLaminar - nIonizationFlag; ++ m)
                        {
                            faceFluxes(i, j, k, m) = -oRefReNumber * wallFluxes[m];
#ifdef USE_ErrorDataDebug
                            if (isErrorData(faceFluxes(i, j, k, m)))
                            {
                                int zoneID = grid->GetZoneID();
                                PrintToWindow(" zoneID = ", zoneID, "i = ", i, "j = ", j, "k = ", k, "m = ", m);
                                TK_Exit::ExceptionExit("This is in function CorrectViscousFlux", true);
                            }
#endif
                        }

                        for (int m = 0; m < nTemperatureModel - 1; ++ m)
                        {
                            faceFluxes(i, j, k, nLaminar + nChemical + m) = -oRefReNumber * wallFluxes[nLaminar + nChemical + m];
#ifdef USE_ErrorDataDebug
                            if (isErrorData(faceFluxes(i, j, k, nLaminar + nChemical + m)))
                            {
                                int zoneID = grid->GetZoneID();
                                PrintToWindow(" zoneID = ", zoneID, "i = ", i, "j = ", j, "k = ", k, "nLaminar + nChemical + m = ", nLaminar + nChemical + m);
                                TK_Exit::ExceptionExit("This is in function CorrectViscousFlux", true);
                            }
#endif
                        }
                        //Next grid cell.
                        ++ indexOfCell;
                    }//i end.
                }//j end.
            }//k end.
            ++ indexOfWall; //Next surface region.
        }//nBCRegion end.

        delete[]wallVariables;    wallVariables = NULL;
        delete[]wallFluxes;    wallFluxes = NULL;
        delete[]speciesEnthalpy;    speciesEnthalpy = NULL;
        delete[]densityDiffusivity;    densityDiffusivity = NULL;
        delete[]vibrationEnergy;    vibrationEnergy = NULL;
        delete[]electronEnergy;    electronEnergy = NULL;
    }

    void NSSolverStruct::CompVisFluxnew(Grid *gridIn, Gradient *gradient_uvwt, FieldProxy *fluxProxy, int iSurface)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &graduvwtx = gradient_uvwt->GetGradX()->GetField_STR();
        RDouble4D &graduvwty = gradient_uvwt->GetGradY()->GetField_STR();
        RDouble4D &graduvwtz = gradient_uvwt->GetGradZ()->GetField_STR();

        RDouble4D &flux = fluxProxy->GetField_STR();

        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
        int numberOfSpecies = parameters->GetNumberOfSpecies();

        RDouble refGama = parameters->GetRefGama();

        RDouble oPrandtlLaminar = parameters->GetoPrandtlLaminar();
        RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble oRefReNumber = parameters->GetoRefReNumber();

        int nrokplus = GlobalDataBase::GetIntParaFromDB("nrokplus");

        RDouble4D *q_turb, *aniss;
        if (nrokplus > 0)
        {
            q_turb = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_turb"));
            aniss = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("aniss"));
        }

        RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);

        RDouble *qlr = new RDouble[nEquation];
        RDouble *fslr = new RDouble[numberOfSpecies];
        RDouble *dfsdx = new RDouble[numberOfSpecies];
        RDouble *dfsdy = new RDouble[numberOfSpecies];
        RDouble *dfsdz = new RDouble[numberOfSpecies];
        RDouble *fvis = new RDouble[nEquation];
        RDouble *work = new RDouble[nEquation];

        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        int iFaceStart, iFaceEnd, jFaceStart, jFaceEnd, kFaceStart, kFaceEnd;
        grid->GetFaceIterationIndex(iFaceStart, iFaceEnd, jFaceStart, jFaceEnd, kFaceStart, kFaceEnd, iSurface);

        using namespace IDX;

        for (int k = kFaceStart; k <= kFaceEnd; ++k)
        {
            for (int j = jFaceStart; j <= jFaceEnd; ++j)
            {
                for (int i = iFaceStart; i <= iFaceEnd; ++i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    //! Note that if we use iSurface,there must be no omission.
                    RDouble nx = xfv(il, jl, kl, iSurface);
                    RDouble ny = yfv(il, jl, kl, iSurface);
                    RDouble nz = zfv(il, jl, kl, iSurface);

                    RDouble um = half * (q(il, jl, kl, IU) + q(i, j, k, IU));
                    RDouble vm = half * (q(il, jl, kl, IV) + q(i, j, k, IV));
                    RDouble wm = half * (q(il, jl, kl, IW) + q(i, j, k, IW));

                    RDouble mul = half * (viscousLaminar(il, jl, kl) + viscousLaminar(i, j, k));
                    RDouble mut = half * (viscousTurbulence(il, jl, kl) + viscousTurbulence(i, j, k));
                    //mut = 0.0;
                    RDouble vis = mul + mut;

                    RDouble dudx = graduvwtx(il, jl, kl, 0);
                    RDouble dudy = graduvwty(il, jl, kl, 0);
                    RDouble dudz = graduvwtz(il, jl, kl, 0);

                    RDouble dvdx = graduvwtx(il, jl, kl, 1);
                    RDouble dvdy = graduvwty(il, jl, kl, 1);
                    RDouble dvdz = graduvwtz(il, jl, kl, 1);

                    RDouble dwdx = graduvwtx(il, jl, kl, 2);
                    RDouble dwdy = graduvwty(il, jl, kl, 2);
                    RDouble dwdz = graduvwtz(il, jl, kl, 2);

                    RDouble dtdx = graduvwtx(il, jl, kl, 3);
                    RDouble dtdy = graduvwty(il, jl, kl, 3);
                    RDouble dtdz = graduvwtz(il, jl, kl, 3);

                    RDouble qx = 0.0;
                    RDouble qy = 0.0;
                    RDouble qz = 0.0;

                    RDouble kcp = (mul * oPrandtlLaminar + mut * oPrandtlTurbulence) * specificHeatAtConstantPressure;

                    qx += kcp * dtdx;
                    qy += kcp * dtdy;
                    qz += kcp * dtdz;

                    RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

                    //! Stress components.
                    RDouble txx = vis * (two * dudx - divv2p3);
                    RDouble tyy = vis * (two * dvdy - divv2p3);
                    RDouble tzz = vis * (two * dwdz - divv2p3);
                    RDouble txy = vis * (dudy + dvdx);
                    RDouble txz = vis * (dudz + dwdx);
                    RDouble tyz = vis * (dvdz + dwdy);

                    fvis[IR] = 0.0;
                    fvis[IDX::IRU] = nx * txx + ny * txy + nz * txz;
                    fvis[IDX::IRV] = nx * txy + ny * tyy + nz * tyz;
                    fvis[IDX::IRW] = nx * txz + ny * tyz + nz * tzz;
                    fvis[IDX::IRE] = um * fvis[IDX::IRU] + vm * fvis[IDX::IRV] + wm * fvis[IDX::IRW] + (nx * qx + ny * qy + nz * qz);

                    for (int m = 0; m < nLaminar; ++m)
                    {
                        flux(il, jl, kl, m) = -oRefReNumber * fvis[m];
                    }
                }
            }
        }
        delete []fvis;    fvis = NULL;
        delete []work;    work = NULL;
        delete []qlr;    qlr = NULL;
        delete []fslr;    fslr = NULL;
        delete []dfsdx;    dfsdx = NULL;
        delete []dfsdy;    dfsdy = NULL;
        delete []dfsdz;    dfsdz = NULL;
    }

    void NSSolverStruct::CompVisFluxnew(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &gradqx = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
        RDouble4D &gradqy = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
        RDouble4D &gradqz = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

        RDouble4D &gradtx = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceX"));
        RDouble4D &gradty = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceY"));
        RDouble4D &gradtz = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceZ"));

        RDouble4D &flux = fluxProxy->GetField_STR();
        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
        RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
        RDouble4D &qTurbulence = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_turb"));
        RDouble3D &wallDistance = *grid->GetWallDist();

        RDouble4D &heatTransferCoeff = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("heatTransferCoeff"));
        int wallFunctionType = GlobalDataBase::GetIntParaFromDB("wallFunctionType");

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
        int numberOfSpecies = parameters->GetNumberOfSpecies();
        RDouble refGama = parameters->GetRefGama();
        int nElectronIndex = parameters->GetIndexOfElectron();
        int nTurblenceForChemical = parameters->GetnTurblenceForChemical();
        int nEnergyRecycle = parameters->GetnEnergyRecycle();

        int transitionType = parameters->GetTransitionType();

        RDouble oPrandtlLaminar = parameters->GetoPrandtlLaminar();
        RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble refReNumber = parameters->GetRefReNumber();
        RDouble oRefReNumber = parameters->GetoRefReNumber();

        int nrokplus = GlobalDataBase::GetIntParaFromDB("nrokplus");
        RDouble4D *heatConductivity = 0, *massDiffusionCoef = 0, *speciesEnthalpy = 0, *speciesEvs = 0, *speciesEes = 0;
        RDouble3D *totalCp = 0;
        if (nChemical > 0)
        {
            heatConductivity = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("lambda"));
            massDiffusionCoef = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rhoD"));
            speciesEnthalpy = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("hSpecies"));
            totalCp = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCp"));

            speciesEvs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEvs"));
            speciesEes = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEes"));
        }
        RDouble kcp[3] = { 0 }, heatFlux = 0.0, dcsdx = 0.0, dcsdy = 0.0, dcsdz = 0.0, hs = 0.0, roDs = 0.0;
        RDouble Tv = 0.0, Te = 0.0, energy = 0.0;
        RDouble4D *q_turb = 0, *aniss = 0;
        if (nrokplus > 0)
        {
            q_turb = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_turb"));
            aniss = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("aniss"));
        }

        RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);
        RDouble *qlr = new RDouble[nEquation]();
        RDouble *fslr = new RDouble[numberOfSpecies]();
        RDouble *dfsdx = new RDouble[numberOfSpecies]();
        RDouble *dfsdy = new RDouble[numberOfSpecies]();
        RDouble *dfsdz = new RDouble[numberOfSpecies]();
        RDouble *fvis = new RDouble[nEquation]();
        RDouble *work = new RDouble[nEquation]();
        RDouble *vibrationEnergy = new RDouble[numberOfSpecies]();
        RDouble *electronEnergy = new RDouble[numberOfSpecies]();
        //RDouble primVars[MAX_SPECIES_NUM], speciesEtr[MAX_SPECIES_NUM], speciesH[MAX_SPECIES_NUM], totalEtr, totalEv, totalH0;
        RDouble wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble twRatio = 1.0;
        if (wallTemperature > 0.0)
        {
            twRatio = wallTemperature / (refDimensionalTemperature * (1.0 + 0.2 * refMachNumber * refMachNumber * 0.85));
        }
        twRatio = sqrt(twRatio);

        int iFaceStart = 0, iFaceEnd = 0, jFaceStart = 0, jFaceEnd = 0, kFaceStart = 0, kFaceEnd = 0;
        grid->GetFaceIterationIndex(iFaceStart, iFaceEnd, jFaceStart, jFaceEnd, kFaceStart, kFaceEnd, iSurface);

        int il1, jl1, kl1;
        grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

        using namespace IDX;

        for (int k = kFaceStart; k <= kFaceEnd; ++ k)
        {
            for (int j = jFaceStart; j <= jFaceEnd; ++ j)
            {
                for (int i = iFaceStart; i <= iFaceEnd; ++ i)
                {
                    int il, jl, kl;
                    grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, iSurface);

                    //! Note that if we use iSurface,there must be no omission.
                    RDouble nx = xfv(il, jl, kl, iSurface);
                    RDouble ny = yfv(il, jl, kl, iSurface);
                    RDouble nz = zfv(il, jl, kl, iSurface);

                    RDouble um = half * (q(il, jl, kl, IU) + q(i, j, k, IU));
                    RDouble vm = half * (q(il, jl, kl, IV) + q(i, j, k, IV));
                    RDouble wm = half * (q(il, jl, kl, IW) + q(i, j, k, IW));

                    RDouble mul = half * (viscousLaminar(il, jl, kl) + viscousLaminar(i, j, k));
                    RDouble mut = half * (viscousTurbulence(il, jl, kl) + viscousTurbulence(i, j, k));
                    RDouble vis = mul + mut;

                    RDouble dudx = gradqx(il, jl, kl, IU);
                    RDouble dudy = gradqy(il, jl, kl, IU);
                    RDouble dudz = gradqz(il, jl, kl, IU);

                    RDouble dvdx = gradqx(il, jl, kl, IV);
                    RDouble dvdy = gradqy(il, jl, kl, IV);
                    RDouble dvdz = gradqz(il, jl, kl, IV);

                    RDouble dwdx = gradqx(il, jl, kl, IW);
                    RDouble dwdy = gradqy(il, jl, kl, IW);
                    RDouble dwdz = gradqz(il, jl, kl, IW);

                    RDouble dtdx = gradtx(il, jl, kl, ITT);
                    RDouble dtdy = gradty(il, jl, kl, ITT);
                    RDouble dtdz = gradtz(il, jl, kl, ITT);

                    if (nChemical == 0)
                    {
                        kcp[ITT] = (mul * oPrandtlLaminar + mut * oPrandtlTurbulence) * specificHeatAtConstantPressure;

                        //if (wallFunctionType == WALLFUNCTION::STANDARD)
                        //{
                        //    RDouble hfCoeff = heatTransferCoeff(i, j, k, ITT) + heatTransferCoeff(il, jl, kl, ITT);
                        //    if (abs(hfCoeff) > SMALL)
                        //    {
                        //        kcp[ITT] += hfCoeff;
                        //    }
                        //}

                        //if (transitionType == IREGAMA && refMachNumber >= 3.0)
                        //{
                        //    RDouble density =  half * (q(il, jl, kl, IR) + q(i, j, k,IR));
                        //    RDouble qKe   =  half * (qTurbulence(il, jl, kl, IKE)   + qTurbulence(i, j, k,IKE));
                        //    RDouble qGama =  half * (qTurbulence(il, jl, kl, IGAMA) + qTurbulence(i, j, k,IGAMA));
                        //    RDouble qk = qKe * qGama;

                        //    int iLe = MAX(i, 1);
                        //    int jLe = MAX(j, 1);
                        //    int kLe = MAX(k, 1);
                        //    int iRe = MIN(il, ni - 1);
                        //    int jRe = MIN(jl, nj - 1);
                        //    int kRe = MIN(kl, nk - 1);

                        //    RDouble wdist = half * (wallDistance(iLe, jLe, kLe) + wallDistance(iRe, jRe, kRe));
                        //    RDouble var = -sqrt(qk) * wdist * density * refReNumber / mul;
                        //    RDouble vprt = 0.41 / pow(0.3, 1.5) / 3.12 * (1.0 - exp(var / 70.0 / twRatio)) / (1.0 - exp(var / 92.0));
                        //    vprt = MAX(MIN(vprt, 1.0), 0.8);
                        //    kcp[ITT] = (mul * oPrandtlLaminar + mut / vprt) * specificHeatAtConstantPressure;
                        //}
                        heatFlux = kcp[ITT] * (nx * dtdx + ny * dtdy + nz * dtdz);
                    }
                    else
                    {
                        kcp[ITT] = 0.5 * ((*heatConductivity)(i, j, k, ITT) + (*heatConductivity)(il, jl, kl, ITT));
                        if (nTurblenceForChemical == 0)
                        {
                            kcp[ITT] *= (1.0 + mut * oPrandtlTurbulence / (mul * oPrandtlLaminar + SMALL));//! modified by dms.
                        }
                        else
                        {
                            kcp[ITT] += mut * oPrandtlTurbulence * half * ((*totalCp)(il, jl, kl) + (*totalCp)(i, j, k));  //! modified by dms. 
                        }

                        heatFlux = kcp[ITT] * (nx * dtdx + ny * dtdy + nz * dtdz);

                        for (int it = 1; it < nTemperatureModel; ++ it)
                        {
                            dtdx = gradtx(il, jl, kl, it);
                            dtdy = gradty(il, jl, kl, it);
                            dtdz = gradtz(il, jl, kl, it);
                            kcp[it] = 0.5 * ((*heatConductivity)(i, j, k, it) + (*heatConductivity)(il, jl, kl, it)); // * kturb; //! modified by dms. 
                            fvis[nLaminar + it] = kcp[it] * (nx * dtdx + ny * dtdy + nz * dtdz);
                            heatFlux += fvis[nLaminar + it];
                        }

                        //! Obtain the gradient of species.
                        for (int s = 0; s < numberOfSpecies; ++ s)
                        {
                            dcsdx = gradqx(il, jl, kl, nNSEquation + s);
                            dcsdy = gradqy(il, jl, kl, nNSEquation + s);
                            dcsdz = gradqz(il, jl, kl, nNSEquation + s);
                            hs = 0.5 * ((*speciesEnthalpy)(i, j, k, s) + (*speciesEnthalpy)(il, jl, kl, s));
                            roDs = 0.5 * ((*massDiffusionCoef)(i, j, k, s) + (*massDiffusionCoef)(il, jl, kl, s));    //! modified by dms.

                            fvis[nNSEquation + s] = roDs * (nx * dcsdx + ny * dcsdy + nz * dcsdz);
                            heatFlux += hs * fvis[nNSEquation + s];
                        }

                        if (nTemperatureModel == 2)
                        {
                            if (nEnergyRecycle == 0)
                            {
                                Tv = 0.5 * (t(i, j, k, ITV) + t(il, jl, kl, ITV));
                                //! Compute the vibration and electron energies.
                                gas->ComputeVibrationEnergy(Tv, vibrationEnergy);
                                gas->ComputeElectronEnergy(Tv, electronEnergy);

                                for (int s = 0; s < numberOfSpecies; ++ s)
                                {
                                    energy = vibrationEnergy[s] + electronEnergy[s];
                                    fvis[nLaminar + 1] += energy * fvis[nNSEquation + s];
                                }
                            }
                            else
                            {
                                for (int s = 0; s < numberOfSpecies; ++ s)
                                {
                                    energy = (*speciesEvs)(i, j, k, s) + (*speciesEvs)(il, jl, kl, s);
                                    energy += (*speciesEes)(i, j, k, s) + (*speciesEes)(il, jl, kl, s);
                                    fvis[nLaminar + 1] += 0.5 * energy * fvis[nNSEquation + s];
                                }
                            }

                        }
                        else if (nTemperatureModel == 3)
                        {
                            if (nEnergyRecycle == 0)
                            {
                                Tv = 0.5 * (t(i, j, k, ITV) + t(il, jl, kl, ITV));
                                Te = 0.5 * (t(i, j, k, ITE) + t(il, jl, kl, ITE));
                                //! Compute the vibration and electron energies.
                                gas->ComputeVibrationEnergy(Tv, vibrationEnergy);
                                gas->ComputeElectronEnergy(Te, electronEnergy);
                                for (int s = 0; s < numberOfSpecies; ++ s)
                                {
                                    fvis[nLaminar + 1] += vibrationEnergy[s] * fvis[nNSEquation + s];
                                    fvis[nLaminar + 2] += electronEnergy[s] * fvis[nNSEquation + s];
                                }
                            }
                            else
                            {
                                for (int s = 0; s < numberOfSpecies; ++ s)
                                {
                                    fvis[nLaminar + 1] += 0.5 * ((*speciesEvs)(i, j, k, s) + (*speciesEvs)(il, jl, kl, s)) * fvis[nNSEquation + s];
                                    fvis[nLaminar + 2] += 0.5 * ((*speciesEes)(i, j, k, s) + (*speciesEes)(il, jl, kl, s)) * fvis[nNSEquation + s];
                                }
                            }
                        }
                        fvis[nLaminar] = 0.0;
                        if (nElectronIndex >= 0)
                        {
                            fvis[nLaminar - 1] = 0.0;
                        }
                    }

                    RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

                    //! Stress components.
                    RDouble txx = vis * (two * dudx - divv2p3);
                    RDouble tyy = vis * (two * dvdy - divv2p3);
                    RDouble tzz = vis * (two * dwdz - divv2p3);
                    RDouble txy = vis * (dudy + dvdx);
                    RDouble txz = vis * (dudz + dwdx);
                    RDouble tyz = vis * (dvdz + dwdy);

                    fvis[IR] = 0.0;
                    fvis[IDX::IRU] = nx * txx + ny * txy + nz * txz;
                    fvis[IDX::IRV] = nx * txy + ny * tyy + nz * tyz;
                    fvis[IDX::IRW] = nx * txz + ny * tyz + nz * tzz;
                    fvis[IDX::IRE] = um * fvis[IDX::IRU] + vm * fvis[IDX::IRV] + wm * fvis[IDX::IRW] + heatFlux;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        flux(il, jl, kl, m) = -oRefReNumber * fvis[m];
#ifdef USE_ErrorDataDebug
                        if (isErrorData(flux(il, jl, kl, m)))
                        {
                            int zoneID = grid->GetZoneID();
                            PrintToWindow(" zoneID = ", zoneID, "il = ", il, "jl = ", jl, "kl = ", kl, "m = ", m);
                            TK_Exit::ExceptionExit("This is in function CompVisFluxnew", true);
                        }
#endif
                    }
                }
            }
        }

        delete[] fvis;    fvis = NULL;
        delete[] work;    work = NULL;
        delete[] qlr;     qlr = NULL;
        delete[] fslr;    fslr = NULL;
        delete[] dfsdx;    dfsdx = NULL;
        delete[] dfsdy;    dfsdy = NULL;
        delete[] dfsdz;    dfsdz = NULL;
        delete[] vibrationEnergy;    vibrationEnergy = NULL;
        delete[] electronEnergy;    electronEnergy = NULL;
    }

    void NSSolverStruct::CompVisFluxLES(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &gradqx = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
        RDouble4D &gradqy = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
        RDouble4D &gradqz = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

        RDouble4D &gradtx = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceX"));
        RDouble4D &gradty = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceY"));
        RDouble4D &gradtz = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTemperatureFaceZ"));

        RDouble4D &flux = fluxProxy->GetField_STR();
        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));

        RDouble3D &vist = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
        RDouble3D &subgridScaleEnergy = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));
        RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
        int numberOfSpecies = parameters->GetNumberOfSpecies();
        RDouble refGama = parameters->GetRefGama();
        int nElectronIndex = parameters->GetIndexOfElectron();

        RDouble oPrandtlLaminar = parameters->GetoPrandtlLaminar();
        //RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble oRefReNumber = parameters->GetoRefReNumber();

        int nrokplus = GlobalDataBase::GetIntParaFromDB("nrokplus");
        RDouble4D *heatConductivity = 0, *massDiffusionCoef = 0, *speciesEnthalpy = 0;
        if (nChemical > 0)
        {
            heatConductivity = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("lambda"));
            massDiffusionCoef = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rhoD"));
            speciesEnthalpy = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("hSpecies"));
        }
        RDouble kcp[3] = { 0 }, heatFlux = 0.0, dcsdx = 0.0, dcsdy = 0.0, dcsdz = 0.0, hs = 0.0, roDs = 0.0;
        RDouble Tv = 0.0, Te = 0.0, energy = 0.0;
        RDouble4D *q_turb = 0, *aniss = 0;
        if (nrokplus > 0)
        {
            q_turb = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_turb"));
            aniss = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("aniss"));
        }

        RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);
        RDouble *qlr = new RDouble[nEquation]();
        RDouble *fslr = new RDouble[numberOfSpecies]();
        RDouble *dfsdx = new RDouble[numberOfSpecies]();
        RDouble *dfsdy = new RDouble[numberOfSpecies]();
        RDouble *dfsdz = new RDouble[numberOfSpecies]();
        RDouble *fvis = new RDouble[nEquation]();
        RDouble *work = new RDouble[nEquation]();
        RDouble *vibrationEnergy = new RDouble[numberOfSpecies]();
        RDouble *electronEnergy = new RDouble[numberOfSpecies]();

        int iFaceStart = 0, iFaceEnd = 0, jFaceStart = 0, jFaceEnd = 0, kFaceStart = 0, kFaceEnd = 0;
        grid->GetFaceIterationIndex(iFaceStart, iFaceEnd, jFaceStart, jFaceEnd, kFaceStart, kFaceEnd, iSurface);

        int il1, jl1, kl1;
        grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

        using namespace IDX;

        string sgsmodel = " ";
        GlobalDataBase::GetData("sgsmodel", &sgsmodel, PHSTRING, 1);

        int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

        RDouble viscousTurbulenceMaximum = -1.0e30;
        RDouble viscousTurbulenceMinimum = 1.0e30;

        int imax = 1;
        int jmax = 1;
        int kmax = 1;

        int imin = 1;
        int jmin = 1;
        int kmin = 1;

        for (int k = kFaceStart; k <= kFaceEnd; ++ k)
        {
            for (int j = jFaceStart; j <= jFaceEnd; ++ j)
            {
                for (int i = iFaceStart; i <= iFaceEnd; ++ i)
                {
                    int il, jl, kl;
                    grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, iSurface);

                    //! Note that if we use iSurface,there must be no omission.
                    RDouble nx = xfv(il, jl, kl, iSurface);
                    RDouble ny = yfv(il, jl, kl, iSurface);
                    RDouble nz = zfv(il, jl, kl, iSurface);

                    RDouble um = half * (q(il, jl, kl, IU) + q(i, j, k, IU));
                    RDouble vm = half * (q(il, jl, kl, IV) + q(i, j, k, IV));
                    RDouble wm = half * (q(il, jl, kl, IW) + q(i, j, k, IW));

                    RDouble mul = half * (viscousLaminar(il, jl, kl) + viscousLaminar(i, j, k));

                    RDouble dudx = gradqx(il, jl, kl, IU);
                    RDouble dudy = gradqy(il, jl, kl, IU);
                    RDouble dudz = gradqz(il, jl, kl, IU);

                    RDouble dvdx = gradqx(il, jl, kl, IV);
                    RDouble dvdy = gradqy(il, jl, kl, IV);
                    RDouble dvdz = gradqz(il, jl, kl, IV);

                    RDouble dwdx = gradqx(il, jl, kl, IW);
                    RDouble dwdy = gradqy(il, jl, kl, IW);
                    RDouble dwdz = gradqz(il, jl, kl, IW);

                    RDouble dtdx = gradtx(il, jl, kl, ITT);
                    RDouble dtdy = gradty(il, jl, kl, ITT);
                    RDouble dtdz = gradtz(il, jl, kl, ITT);

                    RDouble oPrandtlTurbulence = 2.0 / (turbulentPrandtlNumber(il, jl, kl) + turbulentPrandtlNumber(i, j, k));

                    RDouble mut = half * (vist(il, jl, kl) + vist(i, j, k));

                    RDouble hfxSGS = -mut * specificHeatAtConstantPressure * oPrandtlTurbulence * dtdx;
                    RDouble hfySGS = -mut * specificHeatAtConstantPressure * oPrandtlTurbulence * dtdy;
                    RDouble hfzSGS = -mut * specificHeatAtConstantPressure * oPrandtlTurbulence * dtdz;

                    if (nChemical == 0)
                    {
                        kcp[ITT] = (mul * oPrandtlLaminar) * specificHeatAtConstantPressure;
                        heatFlux = kcp[ITT] * (nx * dtdx + ny * dtdy + nz * dtdz) - (nx * hfxSGS + ny * hfySGS + nz * hfzSGS);
                    }
                    else
                    {
                        kcp[ITT] = 0.5 * ((*heatConductivity)(i, j, k, ITT) + (*heatConductivity)(il, jl, kl, ITT));
                        heatFlux = kcp[ITT] * (nx * dtdx + ny * dtdy + nz * dtdz);
                        for (int it = 1; it < nTemperatureModel; ++ it)
                        {
                            dtdx = gradtx(il, jl, kl, it);
                            dtdy = gradty(il, jl, kl, it);
                            dtdz = gradtz(il, jl, kl, it);
                            kcp[it] = 0.5 * ((*heatConductivity)(i, j, k, it) + (*heatConductivity)(il, jl, kl, it));
                            fvis[nLaminar + it] = kcp[it] * (nx * dtdx + ny * dtdy + nz * dtdz);
                            heatFlux += fvis[nLaminar + it];
                        }

                        //! Obtain the gradient of species.
                        for (int s = 0; s < numberOfSpecies; ++ s)
                        {
                            dcsdx = gradqx(il, jl, kl, nNSEquation + s);
                            dcsdy = gradqy(il, jl, kl, nNSEquation + s);
                            dcsdz = gradqz(il, jl, kl, nNSEquation + s);
                            hs = 0.5 * ((*speciesEnthalpy)(i, j, k, s) + (*speciesEnthalpy)(il, jl, kl, s));
                            roDs = 0.5 * ((*massDiffusionCoef)(i, j, k, s) + (*massDiffusionCoef)(il, jl, kl, s));

                            fvis[nNSEquation + s] = roDs * (nx * dcsdx + ny * dcsdy + nz * dcsdz);
                            heatFlux += hs * fvis[nNSEquation + s];
                        }

                        if (nTemperatureModel == 2)
                        {
                            Tv = 0.5 * (t(i, j, k, ITV) + t(il, jl, kl, ITV));
                            //! Compute the vibration and electron energies.
                            gas->ComputeVibrationEnergy(Tv, vibrationEnergy);
                            gas->ComputeElectronEnergy(Tv, electronEnergy);
                            for (int s = 0; s < numberOfSpecies; ++ s)
                            {
                                energy = vibrationEnergy[s] + electronEnergy[s];
                                fvis[nLaminar + 1] += energy * fvis[nNSEquation + s];
                            }
                        }
                        else if (nTemperatureModel == 3)
                        {
                            Tv = 0.5 * (t(i, j, k, ITV) + t(il, jl, kl, ITV));
                            Te = 0.5 * (t(i, j, k, ITE) + t(il, jl, kl, ITE));
                            //! Compute the vibration and electron energies.
                            gas->ComputeVibrationEnergy(Tv, vibrationEnergy);
                            gas->ComputeElectronEnergy(Te, electronEnergy);
                            for (int s = 0; s < numberOfSpecies; ++ s)
                            {
                                fvis[nLaminar + 1] += vibrationEnergy[s] * fvis[nNSEquation + s];
                                fvis[nLaminar + 2] += electronEnergy[s] * fvis[nNSEquation + s];
                            }
                        }
                        fvis[nLaminar] = 0.0;
                        if (nElectronIndex >= 0)
                        {
                            fvis[nLaminar - 1] = 0.0;
                        }
                    }

                    RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

                    //! viscous stress and subgrid scale stress
                    RDouble txx = mul * (two * dudx - divv2p3);
                    RDouble tyy = mul * (two * dvdy - divv2p3);
                    RDouble tzz = mul * (two * dwdz - divv2p3);
                    RDouble txy = mul * (dudy + dvdx);
                    RDouble txz = mul * (dudz + dwdx);
                    RDouble tyz = mul * (dvdz + dwdy);

                    //RDouble ci = 0.5 * (isotropicConstant(il, jl, kl) + isotropicConstant(i, j, k));
                    //RDouble kSGS = 2.0 * rhom * ci * ss * ss * refReNumber;
                    RDouble kSGS = half * (subgridScaleEnergy(il, jl, kl) + subgridScaleEnergy(i, j, k));

                    RDouble txxSGS = -mut * (two * dudx - divv2p3) + 1.0 / 3.0 * kSGS;
                    RDouble tyySGS = -mut * (two * dvdy - divv2p3) + 1.0 / 3.0 * kSGS;
                    RDouble tzzSGS = -mut * (two * dwdz - divv2p3) + 1.0 / 3.0 * kSGS;
                    RDouble txySGS = -mut * (dudy + dvdx);
                    RDouble txzSGS = -mut * (dudz + dwdx);
                    RDouble tyzSGS = -mut * (dvdz + dwdy);

                    fvis[IR] = 0.0;
                    fvis[IDX::IRU] = nx * (txx - txxSGS) + ny * (txy - txySGS) + nz * (txz - txzSGS);
                    fvis[IDX::IRV] = nx * (txy - txySGS) + ny * (tyy - tyySGS) + nz * (tyz - tyzSGS);
                    fvis[IDX::IRW] = nx * (txz - txzSGS) + ny * (tyz - tyzSGS) + nz * (tzz - tzzSGS);
                    fvis[IDX::IRE] = um * fvis[IRU] + vm * fvis[IRV] + wm * fvis[IRW] + heatFlux;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        flux(il, jl, kl, m) = -oRefReNumber * fvis[m];
                    }

                    if (viscousTurbulenceMaximum < mut)
                    {
                        viscousTurbulenceMaximum = mut;
                        imax = i;
                        jmax = j;
                        kmax = k;
                    }

                    if (viscousTurbulenceMinimum > mut)
                    {
                        viscousTurbulenceMinimum = mut;
                        imin = i;
                        jmin = j;
                        kmin = k;
                    }

                }
            }
        }

        RDouble vist_max = viscousTurbulenceMaximum;
        RDouble vist_min = viscousTurbulenceMinimum;

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
                cout << "vistmax = " << viscousTurbulenceMaximum << " " << imax << " " << jmax << " " << kmax << "\n";
                cout << "vistmin = " << viscousTurbulenceMinimum << " " << imin << " " << jmin << " " << kmin << "\n";
            }
#ifdef PH_PARALLEL
        }
#endif

        delete[]fvis;    fvis = NULL;
        delete[]work;    work = NULL;
        delete[]qlr;    qlr = NULL;
        delete[]fslr;    fslr = NULL;
        delete[]dfsdx;    dfsdx = NULL;
        delete[]dfsdy;    dfsdy = NULL;
        delete[]dfsdz;    dfsdz = NULL;
        delete[]vibrationEnergy;    vibrationEnergy = NULL;
        delete[]electronEnergy;    electronEnergy = NULL;
    }

    void NSSolverStruct::LoadFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &flux = fluxProxy->GetField_STR();

        int nEquation = GetNumberOfEquations();

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));

        int il1, jl1, kl1;
        grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

        Range I, J, K;
        GetRange(ni, nj, nk, 0, -1, I, J, K);

        Range M(0, nEquation - 1);

        residual(I, J, K, M) -= (flux(I + il1, J + jl1, K + kl1, M) - flux(I, J, K, M));
    }

    void NSSolverStruct::UpdateQlQrOnlyforMixGrid(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        FieldProxy *qlProxy = new FieldProxy();
        FieldProxy *qrProxy = new FieldProxy();

        int nDim = GetDim();
        for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
        {
            if (iSurface == 1)
            {
                qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1")));
                qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1")));
            }
            else if (iSurface == 2)
            {
                qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2")));
                qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2")));
            }
            else
            {
                qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3")));
                qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3")));
            }

            //! Apply the limiter to compute left-side variable QL and right-side variable QR.
            GetInvFaceVar(grid, qlProxy, qrProxy, iSurface);

            //! Modify the variables of the wall boundary faces-clz begin.
            CorrectFaceVar(grid, qlProxy, qrProxy, iSurface);
        }

        delete qlProxy;
        delete qrProxy;
    }

    void NSSolverStruct::UpdateInviscidfluxOnlyforMixGrid(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        int structureScheme = GlobalDataBase::GetIntParaFromDB("str_scheme");
        INVScheme inviscidScheme = GetINVScheme(structureScheme);

        FieldProxy *qlProxy = new FieldProxy();
        FieldProxy *qrProxy = new FieldProxy();
        FieldProxy *fluxProxy = new FieldProxy();

        int nDim = GetDim();
        for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
        {
            if (iSurface == 1)
            {
                qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1")));
                qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1")));
                fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxI")));
            }
            else if (iSurface == 2)
            {
                qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2")));
                qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2")));
                fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxJ")));
            }
            else
            {
                qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3")));
                qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3")));
                fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxK")));
            }

            //! To compute the inviscid flux of the faces using the QL and QR.
            InviscidFluxWrap(grid, qlProxy, qrProxy, fluxProxy, inviscidScheme, iSurface);
        }

        delete qlProxy;
        delete qrProxy;
        delete fluxProxy;
    }

    void NSSolverStruct::UpdateViscousfluxOnlyforMixGrid(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int viscousType = parameters->GetViscousType();
        if (viscousType == INVISCID) return;

        FieldProxy *fluxProxy = new FieldProxy();

        int nDim = GetDim();
        for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
        {
            if (iSurface == 1)
            {
                fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxI")));
            }
            else if (iSurface == 2)
            {
                fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxJ")));
            }
            else
            {
                fluxProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxK")));
            }

            ReconGradVolWeightedwithCorrection(gridIn, iSurface);

            CompVisFluxnew(grid, fluxProxy, iSurface);
        }
        delete fluxProxy;
    }

    void NSSolverStruct::LoadGlobalfacefluxtoResOnlyforMixGrid(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));

        RDouble4D &faceInviscidfluxI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxI"));
        RDouble4D &faceInviscidfluxJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxJ"));
        RDouble4D &faceInviscidfluxK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxK"));

        int ist = 1;
        int jst = 1;
        int kst = 1;
        int ied = ni - 1;
        int jed = nj - 1;
        int ked = nk - 1;

        if (nk == 1)
        {
            ked = 1;
        }

        for (int m = 0; m < 5; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        residual(i, j, k, m) -= faceInviscidfluxI(i + 1, j, k, m) - faceInviscidfluxI(i, j, k, m);
                        residual(i, j, k, m) -= faceInviscidfluxJ(i, j + 1, k, m) - faceInviscidfluxJ(i, j, k, m);
                        if (nk == 1) continue;
                        residual(i, j, k, m) -= faceInviscidfluxK(i, j, k + 1, m) - faceInviscidfluxK(i, j, k, m);
                    }
                }
            }
        }

        Param_NSSolverStruct *parameters = GetControlParameters();
        int viscousType = parameters->GetViscousType();
        if (viscousType == INVISCID) return;

        RDouble4D &faceViscousfluxI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxI"));
        RDouble4D &faceViscousfluxJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxJ"));
        RDouble4D &faceViscousfluxK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxK"));

        for (int m = 0; m < 5; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        residual(i, j, k, m) -= faceViscousfluxI(i + 1, j, k, m) - faceViscousfluxI(i, j, k, m);
                        residual(i, j, k, m) -= faceViscousfluxJ(i, j + 1, k, m) - faceViscousfluxJ(i, j, k, m);
                        if (nk == 1) continue;
                        residual(i, j, k, m) -= faceViscousfluxK(i, j, k + 1, m) - faceViscousfluxK(i, j, k, m);
                    }
                }
            }
        }
    }

    void NSSolverStruct::CorrectFaceVar(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface)
    {
        Param_NSSolverStruct *parameters = GetControlParameters();
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &ql = qlProxy->GetField_STR();
        RDouble4D &qr = qrProxy->GetField_STR();
        RDouble4D &primitveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
        RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
        RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());

        RDouble4D &faceVelocityX = *(grid->GetFaceVelocityX());
        RDouble4D &faceVelocityY = *(grid->GetFaceVelocityY());
        RDouble4D &faceVelocityZ = *(grid->GetFaceVelocityZ());

        int viscousType = parameters->GetViscousType();
        RDouble wallTemperature = parameters->GetWallTemperature();
        RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
        RDouble temperatureWall = wallTemperature / refDimensionalTemperature;

        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        RDouble3D *surfaceSlipVariables = NULL;
        if (nSlipBCModel > 0)
        {
            surfaceSlipVariables = reinterpret_cast <RDouble3D *> (gridIn->GetDataPtr("surfaceSlipVariables"));
        }

        using namespace GAS_SPACE;

        RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();
        int nEquation = GetNumberOfEquations();
        RDouble *wallVariables = new RDouble[nEquation];

        using namespace IDX;

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        int nWallBC = 0, nIndexOfCell = 0;

        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int nBCType = structBC->GetBCType();

            if (!IsWall(nBCType))
            {
                continue;
            }

            RDouble uWall = 0.0;
            RDouble vWall = 0.0;
            RDouble wWall = 0.0;

            vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();

            if (globalBCList != 0)
            {
                vector<SimpleBC *>::iterator iter;
                for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
                {
                    SimpleBC *boundaryCondition = *iter;
                    string BCName = structBC->GetBCName();

                    if (BCName != boundaryCondition->GetBCName())
                    {
                        continue;
                    }

                    if (boundaryCondition->CheckParamData("uWall"))
                    {
                        boundaryCondition->GetParamData("uWall", &uWall, PHDOUBLE, 1);
                    }

                    if (boundaryCondition->CheckParamData("vWall"))
                    {
                        boundaryCondition->GetParamData("vWall", &vWall, PHDOUBLE, 1);
                    }

                    if (boundaryCondition->CheckParamData("wWall"))
                    {
                        boundaryCondition->GetParamData("wWall", &wWall, PHDOUBLE, 1);
                    }
                }
            }

            int *faceDirectionIndex = structBC->GetFaceDirectionIndex();
            int directionIndex = structBC->GetFaceDirection() + 1;

            int ist, ied, jst, jed, kst, ked;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            if (directionIndex != iSurface)
            {
                continue;
            }

            //! The following treatment places subscripts on BCFace.
            //! Modified by clz : 2012-7-10.
            int id, jd, kd;
            GetBCFaceIDX(faceDirectionIndex, id, jd, kd);
            kst = kst + kd;
            ked = ked + kd;
            jst = jst + jd;
            jed = jed + jd;
            ist = ist + id;
            ied = ied + id;

            int di, dj, dk;
            di = ABS(faceDirectionIndex[0]);
            dj = ABS(faceDirectionIndex[1]);
            dk = ABS(faceDirectionIndex[2]);
            nIndexOfCell = 0;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        int il = i - di;
                        int jl = j - dj;
                        int kl = k - dk;
                        if (nSlipBCModel > 0)
                        {
                            uWall = (*surfaceSlipVariables)(nWallBC, nIndexOfCell, 3);
                            vWall = (*surfaceSlipVariables)(nWallBC, nIndexOfCell, 4);
                            wWall = (*surfaceSlipVariables)(nWallBC, nIndexOfCell, 5);
                        }

                        RDouble velocityXWall = uWall;
                        RDouble velocityYWall = vWall;
                        RDouble velocityZWall = wWall;
                        if (isUnsteady && isAle)
                        {
                            velocityXWall += faceVelocityX(i, j, k, iSurface);
                            velocityYWall += faceVelocityY(i, j, k, iSurface);
                            velocityZWall += faceVelocityZ(i, j, k, iSurface);
                        }

                        if (viscousType == INVISCID)
                        {
                            RDouble nxs = faceNormalComponentX(i, j, k, iSurface);
                            RDouble nys = faceNormalComponentY(i, j, k, iSurface);
                            RDouble nzs = faceNormalComponentZ(i, j, k, iSurface);

                            RDouble vx = primitveVariables(i - id, j - jd, k - kd, IU);
                            RDouble vy = primitveVariables(i - id, j - jd, k - kd, IV);
                            RDouble vz = primitveVariables(i - id, j - jd, k - kd, IW);
                            RDouble vn = nxs * vx + nys * vy + nzs * vz;

                            ql(il, jl, kl, IR) = primitveVariables(i - id, j - jd, k - kd, IR);
                            ql(il, jl, kl, IU) = primitveVariables(i - id, j - jd, k - kd, IU) - nxs * vn;
                            ql(il, jl, kl, IV) = primitveVariables(i - id, j - jd, k - kd, IV) - nys * vn;
                            ql(il, jl, kl, IW) = primitveVariables(i - id, j - jd, k - kd, IW) - nzs * vn;
                            ql(il, jl, kl, IP) = primitveVariables(i - id, j - jd, k - kd, IP);

                            for (int m = 0; m < nEquation; ++ m)
                            {
                                qr(i, j, k, m) = ql(il, jl, kl, m);
                            }
                        }
                        else if (wallTemperature <= 0.0)
                        {
                            //! Viscous wall, adiabatic.
                            ql(il, jl, kl, IR) = primitveVariables(i - id, j - jd, k - kd, IR);
                            ql(il, jl, kl, IU) = velocityXWall;
                            ql(il, jl, kl, IV) = velocityYWall;
                            ql(il, jl, kl, IW) = velocityZWall;
                            ql(il, jl, kl, IP) = primitveVariables(i - id, j - jd, k - kd, IP);

                            for (int m = 0; m < nEquation; ++ m)
                            {
                                qr(i, j, k, m) = ql(il, jl, kl, m);
                            }
                        }
                        else
                        {
                            //! Isothermal wall.
                            temperatureWall = wallTemperature / refDimensionalTemperature;
                            wallVariables[IP] = primitveVariables(i - id, j - jd, k - kd, IP);
                            wallVariables[IU] = velocityXWall;
                            wallVariables[IV] = velocityYWall;
                            wallVariables[IW] = velocityZWall;

                            RDouble massReciprocal = one;
                            RDouble presureWall = wallVariables[IP];
                            RDouble densityWall = presureWall / (coefficientofstateEquation * temperatureWall * massReciprocal);
                            wallVariables[IR] = densityWall;
                            if (nSlipBCModel > 0)
                            {
                                temperatureWall = (*surfaceSlipVariables)(nWallBC, nIndexOfCell, ITT);
                                wallVariables[IR] = presureWall / (coefficientofstateEquation * temperatureWall * massReciprocal);
                            }

                            for (int m = 0; m < nEquation; ++ m)
                            {
                                ql(il, jl, kl, m) = wallVariables[m];
                                qr(i, j, k, m) = wallVariables[m];
                            }
                        }
                        ++ nIndexOfCell; //Next grid cell.
                    }
                }
            }
            ++ nWallBC; //Next surface.
        }

        delete[]wallVariables;
    }
    //clz end

    void NSSolverStruct::CorrectChemicalFaceVar(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface)
    {
        Param_NSSolverStruct *parameters = GetControlParameters();
        int nChemical = parameters->GetChemicalFlag();

        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &ql = qlProxy->GetField_STR();
        RDouble4D &qr = qrProxy->GetField_STR();
        RDouble4D &primitveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
        RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
        RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());

        RDouble4D &faceVelocityX = *(grid->GetFaceVelocityX());
        RDouble4D &faceVelocityY = *(grid->GetFaceVelocityY());
        RDouble4D &faceVelocityZ = *(grid->GetFaceVelocityZ());

        int viscousType = parameters->GetViscousType();
        RDouble wallTemperature = parameters->GetWallTemperature();
        RDouble2D *surfaceTemperature = NULL;
        RDouble2D *injectionVelocity = NULL;
        if (wallTemperature >= 0.0)
        {
            surfaceTemperature = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("surfaceTemperature"));
        }

        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        RDouble3D *surfaceSlipVariables = NULL;
        if (nSlipBCModel > 0)
        {
            surfaceSlipVariables = reinterpret_cast <RDouble3D *> (gridIn->GetDataPtr("surfaceSlipVariables"));
        }

        int isInjection = GlobalDataBase::GetIntParaFromDB("isInjection");
        if (isInjection == 1)
        {
            injectionVelocity = reinterpret_cast <RDouble2D *> (grid->GetDataPtr("injectionVelocity"));
        }

        RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
        RDouble catalyticCoef = parameters->GetCatalyticCoef();//! The type of catalytic wall condition.

        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

        using namespace GAS_SPACE;

        RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();
        int nSpeciesNumber = gas->GetNumberOfSpecies();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = GetNumberOfEquations();
        int nEnergyRecycle = parameters->GetnEnergyRecycle();

        RDouble *wallVariables = new RDouble[nEquation]();
        RDouble3D *surfaceMassFraction = 0;
        if (nChemical > 0 && nSpeciesNumber > 0)    //! The fully catalytic wall condition.
        {
            surfaceMassFraction = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
        }
        RDouble massReciprocal = one;
        RDouble twall = wallTemperature / refDimensionalTemperature;
        RDouble temperatureWall[3] = { twall, twall, twall };

        using namespace PHENGLEI;
        using namespace IDX;

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        int indexOfWall = 0, indexOfCell = 0;

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int nBCType = structBC->GetBCType();

            if (!IsWall(nBCType) && nBCType != PHENGLEI::ABLATION_SURFACE)
            {
                continue;
            }
            RDouble uWall = 0.0;
            RDouble vWall = 0.0;
            RDouble wWall = 0.0;

            int nSurface = structBC->GetFaceDirection() + 1;

            vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();

            if (globalBCList != 0)
            {
                vector<SimpleBC *>::iterator iter;
                for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
                {
                    SimpleBC *boundaryCondition = *iter;
                    string BCName = structBC->GetBCName();

                    if (BCName != boundaryCondition->GetBCName())
                    {
                        continue;
                    }

                    if (boundaryCondition->CheckParamData("uWall"))
                    {
                        boundaryCondition->GetParamData("uWall", &uWall, PHDOUBLE, 1);
                    }

                    if (boundaryCondition->CheckParamData("vWall"))
                    {
                        boundaryCondition->GetParamData("vWall", &vWall, PHDOUBLE, 1);
                    }

                    if (boundaryCondition->CheckParamData("wWall"))
                    {
                        boundaryCondition->GetParamData("wWall", &wWall, PHDOUBLE, 1);
                    }
                }
            }

            int *faceDirectionIndex = structBC->GetFaceDirectionIndex();
            int directionIndex = structBC->GetFaceDirection() + 1;
            int leftOrRightIndex = structBC->GetFaceLeftOrRightIndex();

            int ist, ied, jst, jed, kst, ked;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            if (directionIndex != iSurface)
            {
                continue;
            }

            //! The following treatment places subscripts on BCFace.
            //! Modified by clz : 2012-7-10.
            int id, jd, kd;
            GetBCFaceIDX(faceDirectionIndex, id, jd, kd);
            kst = kst + kd;
            ked = ked + kd;
            jst = jst + jd;
            jed = jed + jd;
            ist = ist + id;
            ied = ied + id;

            int di, dj, dk;
            di = ABS(faceDirectionIndex[0]);
            dj = ABS(faceDirectionIndex[1]);
            dk = ABS(faceDirectionIndex[2]);

            indexOfCell = 0;
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        int il = i - di;
                        int jl = j - dj;
                        int kl = k - dk;

                        if (isInjection == 1)
                        {
                            int in, jn, kn;
                            structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                            RDouble normalComponetX = -leftOrRightIndex * faceNormalComponentX(in, jn, kn, nSurface);
                            RDouble normalComponetY = -leftOrRightIndex * faceNormalComponentY(in, jn, kn, nSurface);
                            RDouble normalComponetZ = -leftOrRightIndex * faceNormalComponentZ(in, jn, kn, nSurface);

                            RDouble wallVelocity = (*injectionVelocity)(indexOfWall, indexOfCell);

                            uWall = wallVelocity * normalComponetX;
                            vWall = wallVelocity * normalComponetY;
                            wWall = wallVelocity * normalComponetZ;
                        }

                        if (nSlipBCModel > 0)
                        {
                            uWall = (*surfaceSlipVariables)(indexOfWall, indexOfCell, 3);
                            vWall = (*surfaceSlipVariables)(indexOfWall, indexOfCell, 4);
                            wWall = (*surfaceSlipVariables)(indexOfWall, indexOfCell, 5);
                        }

                        RDouble velocityXWall = uWall;
                        RDouble velocityYWall = vWall;
                        RDouble velocityZWall = wWall;
                        if (isUnsteady && isAle)
                        {
                            velocityXWall += faceVelocityX(i, j, k, iSurface);
                            velocityYWall += faceVelocityY(i, j, k, iSurface);
                            velocityZWall += faceVelocityZ(i, j, k, iSurface);
                        }

                        if (viscousType == INVISCID)
                        {
                            RDouble nxs = faceNormalComponentX(i, j, k, iSurface);
                            RDouble nys = faceNormalComponentY(i, j, k, iSurface);
                            RDouble nzs = faceNormalComponentZ(i, j, k, iSurface);

                            RDouble vx = primitveVariables(i - id, j - jd, k - kd, IU);
                            RDouble vy = primitveVariables(i - id, j - jd, k - kd, IV);
                            RDouble vz = primitveVariables(i - id, j - jd, k - kd, IW);
                            RDouble vn = nxs * vx + nys * vy + nzs * vz;

                            ql(il, jl, kl, IR) = primitveVariables(i - id, j - jd, k - kd, IR);
                            ql(il, jl, kl, IU) = primitveVariables(i - id, j - jd, k - kd, IU) - nxs * vn;
                            ql(il, jl, kl, IV) = primitveVariables(i - id, j - jd, k - kd, IV) - nys * vn;
                            ql(il, jl, kl, IW) = primitveVariables(i - id, j - jd, k - kd, IW) - nzs * vn;
                            ql(il, jl, kl, IP) = primitveVariables(i - id, j - jd, k - kd, IP);
                            for (int s = nNSEquation; s < nEquation; ++ s)
                            {
                                ql(il, jl, kl, s) = primitveVariables(i - id, j - jd, k - kd, s);
                            }

                            for (int m = 0; m < nEquation; ++ m)
                            {
                                qr(i, j, k, m) = ql(il, jl, kl, m);
                            }
                        }
                        else if (wallTemperature < 0.0)
                        {
                            //! Viscous wall, adiabatic.
                            ql(il, jl, kl, IR) = primitveVariables(i - id, j - jd, k - kd, IR);
                            ql(il, jl, kl, IU) = velocityXWall;
                            ql(il, jl, kl, IV) = velocityYWall;
                            ql(il, jl, kl, IW) = velocityZWall;
                            ql(il, jl, kl, IP) = primitveVariables(i - id, j - jd, k - kd, IP);
                            for (int s = nNSEquation; s < nEquation; ++ s)
                            {
                                ql(il, jl, kl, s) = primitveVariables(i - id, j - jd, k - kd, s);
                            }

                            for (int m = 0; m < nEquation; ++ m)
                            {
                                qr(i, j, k, m) = ql(il, jl, kl, m);
                            }
                        }
                        else //! Isothermal wall.
                        {
                            wallVariables[IU] = velocityXWall;
                            wallVariables[IV] = velocityYWall;
                            wallVariables[IW] = velocityZWall;
                            wallVariables[IP] = primitveVariables(i - id, j - jd, k - kd, IP);
                            if (nSlipBCModel > 0)
                            {
                                for (int m = 0; m < 3; ++ m)
                                {
                                    temperatureWall[m] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m);
                                }
                            }
                            else
                            {
                                temperatureWall[ITT] = (*surfaceTemperature)(indexOfWall, indexOfCell);
                                temperatureWall[ITV] = temperatureWall[ITT];
                                temperatureWall[ITE] = temperatureWall[ITT];
                            }

                            if (ABS(catalyticCoef) <= EPSILON && nBCType != PHENGLEI::ABLATION_SURFACE) //Fully non-catalytic wall.
                            {
                                for (int s = nNSEquation; s < nNSEquation + nSpeciesNumber; ++ s)
                                {
                                    wallVariables[s] = primitveVariables(i - id, j - jd, k - kd, s);
                                }
                            }
                            else //Fully catalytic wall or finite catalytic wall.
                            {
                                for (int s = nNSEquation; s < nNSEquation + nSpeciesNumber; ++ s)
                                {
                                    wallVariables[s] = (*surfaceMassFraction)(indexOfWall, indexOfCell, s - nNSEquation);
                                }
                            }

                            if (nSlipBCModel > 0)
                            {
                                for (int m = 0; m < nSpeciesNumber; ++ m)
                                {
                                    wallVariables[nNSEquation + m] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m + 6);
                                }
                            }

                            massReciprocal = gas->ComputeMolecularWeightReciprocal(wallVariables);

                            //! Compute vibration and electron energy.
                            /*if(nEnergyRecycle /= 2)
                            {*/
                            if (nTemperatureModel == 2)
                            {
                                wallVariables[nLaminar + 1] = gas->GetMixedGasVibrationEnergy(wallVariables, temperatureWall[ITV]);
                                wallVariables[nLaminar + 1] += gas->GetMixedGasElectronEnergy(wallVariables, temperatureWall[ITV]);
                            }
                            else if (nTemperatureModel == 3)
                            {
                                wallVariables[nLaminar + 1] = gas->GetMixedGasVibrationEnergy(wallVariables, temperatureWall[ITV]);
                                wallVariables[nLaminar + 2] = gas->GetMixedGasElectronEnergy(wallVariables, temperatureWall[ITE]);
                            }
                            /*}
                            else
                            {
                                 if (nTemperatureModel == 2)
                                {
                                    wallVariables[nLaminar + 1]  = temperatureWall[ITV];
                                }
                                else if (nTemperatureModel == 3)
                                {
                                    wallVariables[nLaminar + 1] = temperatureWall[ITV];
                                    wallVariables[nLaminar + 2] = temperatureWall[ITE];
                                }

                            }*/

                            wallVariables[IR] = wallVariables[IP] / (coefficientofstateEquation * temperatureWall[ITT] * massReciprocal);

                            //! Modify the face variables.
                            for (int m = 0; m < nEquation; ++ m)
                            {
                                ql(il, jl, kl, m) = wallVariables[m];
                                qr(i, j, k, m) = wallVariables[m];
                            }

                        }
                        //Next grid cell.
                        ++ indexOfCell;
                    }//i end
                }//j end
            }//k end
            //Next surface region.
            ++ indexOfWall;
        }//iBCRegion end

        delete[] wallVariables;
    }

    void NSSolverStruct::CalPressureFactor(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("rtem"));

        RDouble pressureSensorCoefficient = third;  // for 3-D

        int istp = 1;
        int jstp = 1;
        int kstp = 1;

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);
        int nk = grid->GetNK();

        //! To degenerate for two dimensional problems, modified by eric 20120209.
        if (nk == 1)
        {
            kstp = 0;
            kCellStart = 1;
            kCellEnd = 1;
            pressureSensorCoefficient = half;
        }

        using namespace IDX;
        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)    //! Corrected by clz, 2011-12-20.
                {
                    int im = i - istp;
                    int jm = j - jstp;
                    int km = k - kstp;

                    int ip = i + istp;
                    int jp = j + jstp;
                    int kp = k + kstp;

                    RDouble pim = q(im, j, k, IP);
                    RDouble pjm = q(i, jm, k, IP);
                    RDouble pkm = q(i, j, km, IP);

                    RDouble pip = q(ip, j, k, IP);
                    RDouble pjp = q(i, jp, k, IP);
                    RDouble pkp = q(i, j, kp, IP);

                    RDouble ppp = q(i, j, k, IP);

                    RDouble dpi = ABS((pip - two * ppp + pim) / (pip + two * ppp + pim));
                    RDouble dpj = ABS((pjp - two * ppp + pjm) / (pjp + two * ppp + pjm));
                    RDouble dpk = ABS((pkp - two * ppp + pkm) / (pkp + two * ppp + pkm));

                    rtem(i, j, k) = pressureSensorCoefficient * (dpi + dpj + dpk);  //
                }
            }
        }
        grid->GhostCell3DExceptInterface(rtem);
        //GhostCell3D(rtem, ni, nj, nk);
    }

    void NSSolverStruct::InviscidFlux(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);

        //CalPressureFactor(grid);

        int nDim = GetDim();

        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Range I, J, K;
        GetRange(ni, nj, nk, -2, 1, I, J, K);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nChemical = parameters->GetChemicalFlag();
        int nEquation = GetNumberOfEquations();
        int nInviscidFluxModify = parameters->GetnInviscidFluxModify();

        Range M(0, nEquation - 1);

        FieldProxy *qlProxy = new FieldProxy();
        FieldProxy *qrProxy = new FieldProxy();


        /*int nEnergyRecycle = parameters->GetnEnergyRecycle();
        int ntmodel = parameters->GetTemperatureModel();
        if(nEnergyRecycle == 2 && ntmodel > 1)
        {
            int nLaminar = parameters->GetLaminarNumber();
            RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
            RDouble4D &q = *reinterpret_cast <RDouble4D *> (gridIn->GetDataPtr("q"));
            RDouble4D &temporaryEve = *reinterpret_cast <RDouble4D *> (gridIn->GetDataPtr("temporaryEve"));
            Range Tmode(1, ntmodel - 1);
            temporaryEve(I, J, K, Tmode) = q(I, J, K, Tmode + nLaminar);
            q(I, J, K, Tmode + nLaminar) = t(I, J, K, Tmode);
        }*/


        for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
        {
            if (iSurface == 1)
            {
                qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1")));
                qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1")));
            }
            else if (iSurface == 2)
            {
                qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2")));
                qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2")));
            }
            else
            {
                qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3")));
                qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3")));
            }

            //! Apply the limiter to compute left-side variable QL and right-side variable QR.
            GetInvFaceVar(grid, qlProxy, qrProxy, iSurface);

            //! Modify the variables of the wall boundary faces-clz begin.
            if (nChemical == 0)
            {
                CorrectFaceVar(grid, qlProxy, qrProxy, iSurface);
            }
            else if (nInviscidFluxModify == 0)
            {
                CorrectChemicalFaceVar(grid, qlProxy, qrProxy, iSurface);
            }
            
            //! To compute the inviscid flux terms of the cells with summation of the fluxes of their own surrounding faces.
            ComputeInviscidFlux(grid, qlProxy, qrProxy, iSurface);
        }

        /*if(nEnergyRecycle == 2  && ntmodel > 1)
        {
            int nLaminar = parameters->GetLaminarNumber();
            RDouble4D &temporaryEve = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("temporaryEve"));
            RDouble4D &q = *reinterpret_cast <RDouble4D *> (gridIn->GetDataPtr("q"));
            Range Tmode(1, ntmodel - 1);
            q(I, J, K, Tmode + nLaminar) = temporaryEve(I, J, K, Tmode);
        }*/

        delete qlProxy;
        delete qrProxy;
    }

    void NSSolverStruct::RegisterCFDSolverInterfaceField()
    {
        NSSolver::RegisterCFDSolverInterfaceField();

        int gridTypeSystem = GlobalDataBase::GetIntParaFromDB("sys_gridtype");
        if (gridTypeSystem == MIXGRID)
        {
            NSSolver::RegisterCFDSolverInterfaceFieldUnstruct();
        }
    }

    void NSSolverStruct::RegisterCFDSolverOversetField()
    {
        NSSolver::RegisterCFDSolverOversetField();
    }

    void NSSolverStruct::ReleaseCFDSolverInterfaceField()
    {
        NSSolver::ReleaseCFDSolverInterfaceField();
    }

    void NSSolverStruct::ZeroResidualOfSpecialCells(Grid *gridIn)  //! Added by Guo Yongheng 20160825.
    {
        Param_NSSolverStruct *parameters = GetControlParameters();
        int isOverset = parameters->GetIsOverLapping();
        if (!isOverset)
        {
            return;
        }

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (gridIn->GetDataPtr("res"));
        StructGrid *grid = PHSPACE::StructGridCast(gridIn);
        int nEquation = GetNumberOfEquations();

        int *iBlank = grid->GetCellTypeContainer();

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        int cellLabel = 0;
        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    if (iBlank[cellLabel] != GENERIC_COLOR)
                    {
                        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
                        {
                            residual(i, j, k, iEquation) = 0.0;
                        }
                    }

                    cellLabel += 1;
                }
            }
        }
        return;
    }

    void NSSolverStruct::GetInvFaceVar(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        Range I, J, K;
        GetRange(ni, nj, nk, -2, 1, I, J, K);

        int nEquation = GetNumberOfEquations();

        Range M(0, nEquation - 1);
        using namespace IDX;

        FieldProxy *qProxy = GetFieldProxy(grid, "q");

        RDouble4D &q = qProxy->GetField_STR();
        RDouble4D &ql = qlProxy->GetField_STR();
        RDouble4D &qr = qrProxy->GetField_STR();

        int il1, jl1, kl1;
        grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

        I = Range(1 - il1, ni - 1 + il1);
        J = Range(1 - jl1, nj - 1 + jl1);
        K = Range(1 - kl1, nk - 1 + kl1);
        if (nk == 1)
        {
            K.setRange(1, 1);
        }

        Param_NSSolverStruct *parameters = GetControlParameters();
        RDouble MUSCLCoefXk = parameters->GetMUSCLCoefXk();
        RDouble MUSCLCoefXb = parameters->GetMUSCLCoefXb();

        int structureLimiter = GlobalDataBase::GetIntParaFromDB("str_limiter");
        int isInIniting = GlobalDataBase::GetIntParaFromDB("isInIniting");

        RDouble c1 = 0.25 * (1.0 - MUSCLCoefXk);
        RDouble c2 = 0.25 * (1.0 + MUSCLCoefXk);

        RDouble4D  dql(I, J, K, M, fortranArray);
        RDouble4D  dqr(I, J, K, M, fortranArray);
        RDouble4D dqlr(I, J, K, M, fortranArray);

        if (!grid->IsFinestGrid())
        {
            //! For coarse grid.
            ql(I, J, K, M) = q(I, J, K, M);
            qr(I, J, K, M) = q(I, J, K, M);

            delete qProxy;
            return;
        }

        if (structureLimiter == ILMT_3rd_SMOOTH)
        {
            dql(I, J, K, M) = q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M);
            dqr(I, J, K, M) = q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M);

            dql(I, J, K, IR) = q(I, J, K, IR) * (dql(I, J, K, IR) / (q(I, J, K, IR) - 0.5 * dql(I, J, K, IR)));
            dqr(I, J, K, IR) = q(I, J, K, IR) * (dqr(I, J, K, IR) / (q(I, J, K, IR) + 0.5 * dqr(I, J, K, IR)));
            dql(I, J, K, IP) = q(I, J, K, IP) * (dql(I, J, K, IP) / (q(I, J, K, IP) - 0.5 * dql(I, J, K, IP)));
            dqr(I, J, K, IP) = q(I, J, K, IP) * (dqr(I, J, K, IP) / (q(I, J, K, IP) + 0.5 * dqr(I, J, K, IP)));

            ql(I, J, K, M) = q(I, J, K, M) + half * Third_Order_SMOOTH(dqr(I, J, K, M), dql(I, J, K, M));
            qr(I, J, K, M) = q(I, J, K, M) - half * Third_Order_SMOOTH(dql(I, J, K, M), dqr(I, J, K, M));
        }
        else if (structureLimiter == ILMT_WENO3_JS)
        {
            //! Define WENO3-JS constants.
            const RDouble epsilon = 1.0e-40;
            const RDouble dk02 = 1.0 / 3.0, dk12 = 2.0 / 3.0;
            RDouble4D beta0(I, J, K, M, fortranArray), beta1(I, J, K, M, fortranArray);
            RDouble4D q02(I, J, K, M, fortranArray), q12(I, J, K, M, fortranArray);
            RDouble4D alpha0(I, J, K, M, fortranArray), alpha1(I, J, K, M, fortranArray);

            //! Smooth coefficient.
            beta0(I, J, K, M) = (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)) * (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M));
            beta1(I, J, K, M) = (q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M)) * (q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M));

            //! Non-linear weights (not normalized)
            alpha0(I, J, K, M) = dk02 / ((beta0(I, J, K, M) + epsilon) * (beta0(I, J, K, M) + epsilon));
            alpha1(I, J, K, M) = dk12 / ((beta1(I, J, K, M) + epsilon) * (beta1(I, J, K, M) + epsilon));

            //! Candidate template
            q02(I, J, K, M) = -1.0 / 2.0 * q(I - il1, J - jl1, K - kl1, M) + 3.0 / 2.0 * q(I, J, K, M);
            q12(I, J, K, M) = 1.0 / 2.0 * q(I, J, K, M) + 1.0 / 2.0 * q(I + il1, J + jl1, K + kl1, M);

            //! Compute ql and Normalize
            ql(I, J, K, M) = (alpha0(I, J, K, M) * q02(I, J, K, M) + alpha1(I, J, K, M) * q12(I, J, K, M)) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));

            //! Smooth coefficient
            beta0(I, J, K, M) = (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M)) * (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M));
            beta1(I, J, K, M) = (q(I - il1, J - jl1, K - kl1, M) - q(I, J, K, M)) * (q(I - il1, J - jl1, K - kl1, M) - q(I, J, K, M));

            //! Non-linear weights (not normalized)
            alpha0(I, J, K, M) = dk02 / ((beta0(I, J, K, M) + epsilon) * (beta0(I, J, K, M) + epsilon));
            alpha1(I, J, K, M) = dk12 / ((beta1(I, J, K, M) + epsilon) * (beta1(I, J, K, M) + epsilon));

            //! Candidate template
            q02(I, J, K, M) = -1.0 / 2.0 * q(I + il1, J + jl1, K + kl1, M) + 3.0 / 2.0 * q(I, J, K, M);
            q12(I, J, K, M) = 1.0 / 2.0 * q(I, J, K, M) + 1.0 / 2.0 * q(I - il1, J - jl1, K - kl1, M);

            //! Compute qr and normalize
            qr(I, J, K, M) = (alpha0(I, J, K, M) * q02(I, J, K, M) + alpha1(I, J, K, M) * q12(I, J, K, M)) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));
        }
        else if (structureLimiter == ILMT_WENN3_PRM211)
        {
            //! Define WENN3_PRM211 constants.
            const RDouble EPSILON = 1.0e-45;
            const RDouble DK02 = 1.0 / 3.0, DK12 = 2.0 / 3.0;
            const RDouble C1 = 1.0, C20L = 7.0e7, C20R = 3.0e6, C21L = 1.0e5, C21R = 3.0e6;
            const int M10 = 5, M11 = 4;

            RDouble4D beta0(I, J, K, M, fortranArray), beta1(I, J, K, M, fortranArray);
            RDouble4D q02(I, J, K, M, fortranArray), q12(I, J, K, M, fortranArray);
            RDouble4D alpha0(I, J, K, M, fortranArray), alpha1(I, J, K, M, fortranArray);
            RDouble4D w0(I, J, K, M, fortranArray), w1(I, J, K, M, fortranArray);

            //! Smooth coefficient
            beta0(I, J, K, M) = 1.0 / 4.0 * (3.0 * q(I, J, K, M) - 4.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M)) * (3.0 * q(I, J, K, M) - 4.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M))
                + 13.0 / 12.0 * (q(I, J, K, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M)) * (q(I, J, K, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M));
            beta1(I, J, K, M) = 1.0 / 4.0 * (3.0 * q(I, J, K, M) - 4.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M)) * (3.0 * q(I, J, K, M) - 4.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M))
                + 13.0 / 12.0 * (q(I, J, K, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M)) * (q(I, J, K, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M));

            //! Non-linear weights (not normalized).
            alpha0(I, J, K, M) = DK02 / ((beta0(I, J, K, M) + EPSILON) * (beta0(I, J, K, M) + EPSILON));
            alpha1(I, J, K, M) = DK12 / ((beta1(I, J, K, M) + EPSILON) * (beta1(I, J, K, M) + EPSILON));

            //! Non-linear weights (not normalized).
            w0(I, J, K, M) = alpha0(I, J, K, M) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));
            w1(I, J, K, M) = alpha1(I, J, K, M) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));

            for (int k = K.first(); k <= K.last(); ++k)
            {
                for (int j = 1 - jl1; j <= nj - 1 + jl1; ++j)
                {
                    for (int i = 1 - il1; i <= ni - 1 + il1; ++i)
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {

                            //! Map function.
                            if (w0(i, j, k, m) <= DK02)
                            {
                                w0(i, j, k, m) = DK02 + (w0(i, j, k, m) - DK02) * (w0(i, j, k, m) - DK02) / ((w0(i, j, k, m) - DK02) + C20L * (w0(i, j, k, m) - DK02) * pow(w0(i, j, k, m), M10) - C1 * w0(i, j, k, m) * w0(i, j, k, m));
                            }
                            else
                            {
                                w0(i, j, k, m) = DK02 + (w0(i, j, k, m) - DK02) * (w0(i, j, k, m) - DK02) / ((w0(i, j, k, m) - DK02) + C20R * (w0(i, j, k, m) - DK02) * pow(1.0 - w0(i, j, k, m), M10) - C1 * (1.0 - w0(i, j, k, m)) * (1.0 - w0(i, j, k, m)));
                            }
                            if (w1(i, j, k, m) <= DK12)
                            {
                                w1(i, j, k, m) = DK12 + (w1(i, j, k, m) - DK12) * (w1(i, j, k, m) - DK12) / ((w1(i, j, k, m) - DK12) + C21L * (w1(i, j, k, m) - DK12) * pow(w1(i, j, k, m), M11) - C1 * w1(i, j, k, m) * w1(i, j, k, m));
                            }
                            else
                            {
                                w1(i, j, k, m) = DK12 + (w1(i, j, k, m) - DK12) * (w1(i, j, k, m) - DK12) / ((w1(i, j, k, m) - DK12) + C21R * (w1(i, j, k, m) - DK12) * pow(1.0 - w1(i, j, k, m), M11) - C1 * (1.0 - w1(i, j, k, m)) * (1.0 - w1(i, j, k, m)));
                            }
                        }
                    }
                }
            }
            //! Candidate template.
            q02(I, J, K, M) = -1.0 / 2.0 * q(I - il1, J - jl1, K - kl1, M) + 3.0 / 2.0 * q(I, J, K, M);
            q12(I, J, K, M) = 1.0 / 2.0 * q(I, J, K, M) + 1.0 / 2.0 * q(I + il1, J + jl1, K + kl1, M);

            //! Compute ql
            ql(I, J, K, M) = (w0(I, J, K, M) * q02(I, J, K, M) + w1(I, J, K, M) * q12(I, J, K, M)) / (w0(I, J, K, M) + w1(I, J, K, M));


            //! Smooth coefficient.
            beta0(I, J, K, M) = 1.0 / 4.0 * (3.0 * q(I, J, K, M) - 4.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M)) * (3.0 * q(I, J, K, M) - 4.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M))
                + 13.0 / 12.0 * (q(I, J, K, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M)) * (q(I, J, K, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M));
            beta1(I, J, K, M) = 1.0 / 4.0 * (3.0 * q(I, J, K, M) - 4.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M)) * (3.0 * q(I, J, K, M) - 4.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M))
                + 13.0 / 12.0 * (q(I, J, K, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M)) * (q(I, J, K, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M));

            //! Non-linear weights(not normalized).
            alpha0(I, J, K, M) = DK02 / ((beta0(I, J, K, M) + EPSILON) * (beta0(I, J, K, M) + EPSILON));
            alpha1(I, J, K, M) = DK12 / ((beta1(I, J, K, M) + EPSILON) * (beta1(I, J, K, M) + EPSILON));

            //! Non-linear weights(normalized).
            w0(I, J, K, M) = alpha0(I, J, K, M) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));
            w1(I, J, K, M) = alpha1(I, J, K, M) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));

            for (int k = K.first(); k <= K.last(); ++k)
            {
                for (int j = 1 - jl1; j <= nj - 1 + jl1; ++j)
                {
                    for (int i = 1 - il1; i <= ni - 1 + il1; ++i)
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {
                            //! Map function.
                            if (w0(i, j, k, m) <= DK02)
                            {
                                w0(i, j, k, m) = DK02 + (w0(i, j, k, m) - DK02) * (w0(i, j, k, m) - DK02) / ((w0(i, j, k, m) - DK02) + C20L * (w0(i, j, k, m) - DK02) * pow(w0(i, j, k, m), M10) - C1 * w0(i, j, k, m) * w0(i, j, k, m));
                            }
                            else
                            {
                                w0(i, j, k, m) = DK02 + (w0(i, j, k, m) - DK02) * (w0(i, j, k, m) - DK02) / ((w0(i, j, k, m) - DK02) + C20R * (w0(i, j, k, m) - DK02) * pow(1.0 - w0(i, j, k, m), M10) - C1 * (1.0 - w0(i, j, k, m)) * (1.0 - w0(i, j, k, m)));
                            }
                            if (w1(i, j, k, m) <= DK12)
                            {
                                w1(i, j, k, m) = DK12 + (w1(i, j, k, m) - DK12) * (w1(i, j, k, m) - DK12) / ((w1(i, j, k, m) - DK12) + C21L * (w1(i, j, k, m) - DK12) * pow(w1(i, j, k, m), M11) - C1 * w1(i, j, k, m) * w1(i, j, k, m));
                            }
                            else
                            {
                                w1(i, j, k, m) = DK12 + (w1(i, j, k, m) - DK12) * (w1(i, j, k, m) - DK12) / ((w1(i, j, k, m) - DK12) + C21R * (w1(i, j, k, m) - DK12) * pow(1.0 - w1(i, j, k, m), M11) - C1 * (1.0 - w1(i, j, k, m)) * (1.0 - w1(i, j, k, m)));
                            }
                        }
                    }
                }
            }

            //! Candidate template.
            q02(I, J, K, M) = -1.0 / 2.0 * q(I + il1, J + jl1, K + kl1, M) + 3.0 / 2.0 * q(I, J, K, M);
            q12(I, J, K, M) = 1.0 / 2.0 * q(I, J, K, M) + 1.0 / 2.0 * q(I - il1, J - jl1, K - kl1, M);
            //! 
            qr(I, J, K, M) = (w0(I, J, K, M) * q02(I, J, K, M) + w1(I, J, K, M) * q12(I, J, K, M)) / (w0(I, J, K, M) + w1(I, J, K, M));
        }
        else if (structureLimiter == ILMT_WENN3_ZM)
        {
            //! Define WENN3_ZM constants.
            const RDouble EPSILON = 1.0e-40;
            const RDouble DK02 = 1.0 / 3.0, DK12 = 2.0 / 3.0;
            const RDouble C1 = 0.1, C2 = 0.1, C30 = 29.0, C31 = 29.0;

            RDouble4D beta0(I, J, K, M, fortranArray), beta1(I, J, K, M, fortranArray), w0(I, J, K, M, fortranArray), w1(I, J, K, M, fortranArray), tau(I, J, K, M, fortranArray);
            RDouble4D q02(I, J, K, M, fortranArray), q12(I, J, K, M, fortranArray), alpha0(I, J, K, M, fortranArray), alpha1(I, J, K, M, fortranArray);;

            //! Smooth coefficient.
            beta0(I, J, K, M) = (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)) * (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M));
            beta1(I, J, K, M) = 1.0 / 4.0 * (3.0 * q(I, J, K, M) - 4.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M)) *
                (3.0 * q(I, J, K, M) - 4.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M))
                + 13.0 / 12.0 * (q(I, J, K, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M)) *
                (q(I, J, K, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M));

            for (int k = K.first(); k <= K.last(); ++k)
            {
                for (int j = 1 - jl1; j <= nj - 1 + jl1; ++j)
                {
                    for (int i = 1 - il1; i <= ni - 1 + il1; ++i)
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {

                            //! Global smooth coefficient.
                            tau(i, j, k, m) = abs(1.0 / 4.0 * (q(i + 2 * il1, j + 2 * jl1, k + 2 * kl1, m) - 3.0 * q(i + il1, j + jl1, k + kl1, m) - 21.0 * q(i, j, k, m) + 23.0 * q(i - il1, j - jl1, k - kl1, m))
                                * (q(i + 2 * il1, j + 2 * jl1, k + 2 * kl1, m) - 3.0 * q(i + il1, j + jl1, k + kl1, m) + 3.0 * q(i, j, k, m) - q(i - il1, j - jl1, k - kl1, m)));

                            //! tau/beta
                            w0(i, j, k, m) = tau(i, j, k, m) / (EPSILON + beta0(i, j, k, m));
                            w1(i, j, k, m) = tau(i, j, k, m) / (EPSILON + beta1(i, j, k, m));

                            //! Map function.
                            if (w0(i, j, k, m) < C30)
                            {
                                w0(i, j, k, m) = w0(i, j, k, m) * w0(i, j, k, m) / (w0(i, j, k, m) + (C2 * w0(i, j, k, m) + C1) * (C30 - w0(i, j, k, m)) * (C30 - w0(i, j, k, m)));
                            }
                            if (w1(i, j, k, m) < C31)
                            {
                                w1(i, j, k, m) = w1(i, j, k, m) * w1(i, j, k, m) / (w1(i, j, k, m) + (C2 * w1(i, j, k, m) + C1) * (C31 - w1(i, j, k, m)) * (C31 - w1(i, j, k, m)));
                            }
                        }
                    }
                }
            }

            //! Non-linear weights (not normalized).
            alpha0(I, J, K, M) = DK02 * (1.0 + w0(I, J, K, M));
            alpha1(I, J, K, M) = DK12 * (1.0 + w1(I, J, K, M));

            //! Candidate template.
            q02(I, J, K, M) = -1.0 / 2.0 * q(I - il1, J - jl1, K - kl1, M) + 3.0 / 2.0 * q(I, J, K, M);
            q12(I, J, K, M) = 1.0 / 2.0 * q(I, J, K, M) + 1.0 / 2.0 * q(I + il1, J + jl1, K + kl1, M);

            //! Compute ql.
            ql(I, J, K, M) = (alpha0(I, J, K, M) * q02(I, J, K, M) + alpha1(I, J, K, M) * q12(I, J, K, M)) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));

            //! Smooth coefficient.
            beta0(I, J, K, M) = (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M)) * (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M));
            beta1(I, J, K, M) = 1.0 / 4.0 * (3.0 * q(I, J, K, M) - 4.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M)) *
                (3.0 * q(I, J, K, M) - 4.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M))
                + 13.0 / 12.0 * (q(I, J, K, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M)) *
                (q(I, J, K, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M));

            for (int k = K.first(); k <= K.last(); ++k)
            {
                for (int j = 1 - jl1; j <= nj - 1 + jl1; ++j)
                {
                    for (int i = 1 - il1; i <= ni - 1 + il1; ++i)
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {

                            //! Global Smooth coefficient.
                            tau(i, j, k, m) = abs(1.0 / 4.0 * (q(i - 2 * il1, j - 2 * jl1, k - 2 * kl1, m) - 3.0 * q(i - il1, j - jl1, k - kl1, m) - 21.0 * q(i, j, k, m) + 23.0 * q(i + il1, j + jl1, k + kl1, m))
                                * (q(i - 2 * il1, j - 2 * jl1, k - 2 * kl1, m) - 3.0 * q(i - il1, j - jl1, k - kl1, m) + 3.0 * q(i, j, k, m) - q(i + il1, j + jl1, k + kl1, m)));

                            //! tau/beta
                            w0(i, j, k, m) = tau(i, j, k, m) / (EPSILON + beta0(i, j, k, m));
                            w1(i, j, k, m) = tau(i, j, k, m) / (EPSILON + beta1(i, j, k, m));

                            //! Map function.
                            if (w0(i, j, k, m) <= C30)
                            {
                                w0(i, j, k, m) = w0(i, j, k, m) * w0(i, j, k, m) / (w0(i, j, k, m) + (C2 * w0(i, j, k, m) + C1) * (C30 - w0(i, j, k, m)) * (C30 - w0(i, j, k, m)));
                            }
                            if (w1(i, j, k, m) <= C31)
                            {
                                w1(i, j, k, m) = w1(i, j, k, m) * w1(i, j, k, m) / (w1(i, j, k, m) + (C2 * w1(i, j, k, m) + C1) * (C31 - w1(i, j, k, m)) * (C31 - w1(i, j, k, m)));
                            }
                        }
                    }
                }
            }

            //! Non-linear weights (not normalized)
            alpha0(I, J, K, M) = DK02 * (1.0 + w0(I, J, K, M));
            alpha1(I, J, K, M) = DK12 * (1.0 + w1(I, J, K, M));

            //! Candidate template.
            q02(I, J, K, M) = -1.0 / 2.0 * q(I + il1, J + jl1, K + kl1, M) + 3.0 / 2.0 * q(I, J, K, M);
            q12(I, J, K, M) = 1.0 / 2.0 * q(I, J, K, M) + 1.0 / 2.0 * q(I - il1, J - jl1, K - kl1, M);

            //! Compute qr.
            qr(I, J, K, M) = (alpha0(I, J, K, M) * q02(I, J, K, M) + alpha1(I, J, K, M) * q12(I, J, K, M)) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));
        }
        else if (structureLimiter == ILMT_WENN3_ZES2)
        {
            //! Define WENN3_ZES2 constants.
            const RDouble EPSILON = 1.0e-40;
            const RDouble DK02 = 1.0 / 3.0, DK12 = 2.0 / 3.0;
            const RDouble cAlpha = 0.15, cBeta1 = 0.15, cBeta0Max = 1.0, cBeta0Min = 1.0e-8, kc = 0.75, psic = 0.3;

            RDouble df1, df2, df3, df4, kn, psi, sigma, kmkc, signkmkc;

            RDouble4D w0(I, J, K, M, fortranArray), w1(I, J, K, M, fortranArray), cbeta0(I, J, K, M, fortranArray);
            RDouble4D q02(I, J, K, M, fortranArray), q12(I, J, K, M, fortranArray), alpha0(I, J, K, M, fortranArray), alpha1(I, J, K, M, fortranArray);
            RDouble4D beta20(I, J, K, M, fortranArray), beta21(I, J, K, M, fortranArray), beta0(I, J, K, M, fortranArray), beta1(I, J, K, M, fortranArray), tau(I, J, K, M, fortranArray);

            beta20(I, J, K, M) = (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)) * (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M));
            beta21(I, J, K, M) = (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M)) * (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M));

            for (int k = K.first(); k <= K.last(); ++k)
            {
                for (int j = 1 - jl1; j <= nj - 1 + jl1; ++j)
                {
                    for (int i = 1 - il1; i <= ni - 1 + il1; ++i)
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {
                            // !Compute cbeta0.
                            df1 = (q(i - 2 * il1, j - 2 * jl1, k - 2 * kl1, m) - 8. * q(i - il1, j - jl1, k - kl1, m) + 8. * q(i + il1, j + jl1, k + kl1, m) - q(i + 2 * il1, j + 2 * jl1, k + 2 * kl1, m)) / 12.;
                            df2 = (-q(i - 2 * il1, j - 2 * jl1, k - 2 * kl1, m) + 16. * q(i - il1, j - jl1, k - kl1, m) - 30. * q(i, j, k, m) + 16. * q(i + il1, j + jl1, k + kl1, m) - q(i + 2 * il1, j + 2 * jl1, k + 2 * kl1, m)) / 12.;
                            df3 = (-q(i - 2 * il1, j - 2 * jl1, k - 2 * kl1, m) + 2. * q(i - il1, j - jl1, k - kl1, m) - 2. * q(i + il1, j + jl1, k + kl1, m) + q(i + 2 * il1, j + 2 * jl1, k + 2 * kl1, m)) / 2.;
                            df4 = (q(i - 2 * il1, j - 2 * jl1, k - 2 * kl1, m) - 4. * q(i - il1, j - jl1, k - kl1, m) + 6. * q(i, j, k, m) - 4. * q(i + il1, j + jl1, k + kl1, m) + q(i + 2 * il1, j + 2 * jl1, k + 2 * kl1, m));
                            kn = sqrt((abs(df3) + abs(df4)) / (abs(df1) + abs(df2) + 1.0e-3));
                            psi = 1. - abs(beta20(i, j, k, m) - beta21(i, j, k, m)) / (beta20(i, j, k, m) + beta21(i, j, k, m) + EPSILON);
                            psi = min(1.0, psi / psic);
                            kmkc = kc - kn;
                            signkmkc = 2.0 * (((signed char*)&kmkc)[sizeof(kmkc) - 1] >> 7) + 1.0;
                            sigma = (1.0 + signkmkc) / 2. + (1.0 - signkmkc) / 2. * psi;
                            cbeta0(i, j, k, m) = (cBeta0Min + sigma * sigma * (cBeta0Max - cBeta0Min));

                        }
                    }
                }
            }

            //! Smooth coefficient.
            beta0(I, J, K, M) = beta20(I, J, K, M) + cbeta0 * (q(I - il1, J - jl1, K - kl1, M) - 2.0 * q(I, J, K, M) + q(I + il1, J + jl1, K + kl1, M)) * (q(I - il1, J - jl1, K - kl1, M) - 2.0 * q(I, J, K, M) + q(I + il1, J + jl1, K + kl1, M));
            beta1(I, J, K, M) = beta21(I, J, K, M) + cBeta1 * (q(I, J, K, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M)) * (q(I, J, K, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M));

            //! Global Smooth coefficient.
            tau(I, J, K, M) = ((2. * q(I + il1, J + jl1, K + kl1, M) - 3. * q(I, J, K, M) + 1. * q(I - il1, J - jl1, K - kl1, M)) * (q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M) - 3. * q(I + il1, J + jl1, K + kl1, M) + 3. * q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)));

            //! tau/beta
            w0(I, J, K, M) = tau / (EPSILON + beta0);
            w1(I, J, K, M) = tau / (EPSILON + beta1);

            //! Non-linear weights(not normalized )
            alpha0(I, J, K, M) = DK02 * (1.0 + cAlpha * w0(I, J, K, M) * w0(I, J, K, M));
            alpha1(I, J, K, M) = DK12 * (1.0 + cAlpha * w1(I, J, K, M) * w1(I, J, K, M));

            //! Candidate template.
            q02(I, J, K, M) = -1.0 / 2.0 * q(I - il1, J - jl1, K - kl1, M) + 3.0 / 2.0 * q(I, J, K, M);
            q12(I, J, K, M) = 1.0 / 2.0 * q(I, J, K, M) + 1.0 / 2.0 * q(I + il1, J + jl1, K + kl1, M);

            //! Compute ql.
            ql(I, J, K, M) = (alpha0(I, J, K, M) * q02(I, J, K, M) + alpha1(I, J, K, M) * q12(I, J, K, M)) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));

            beta20(I, J, K, M) = (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M)) * (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M));
            beta21(I, J, K, M) = (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)) * (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M));

            for (int k = K.first(); k <= K.last(); ++k)
            {
                for (int j = 1 - jl1; j <= nj - 1 + jl1; ++j)
                {
                    for (int i = 1 - il1; i <= ni - 1 + il1; ++i)
                    {
                        for (int m = 0; m < nEquation; ++m)
                        {

                            df1 = (q(i + 2 * il1, j + 2 * jl1, k + 2 * kl1, m) - 8. * q(i + il1, j + jl1, k + kl1, m) + 8. * q(i - il1, j - jl1, k - kl1, m) - q(i - 2 * il1, j - 2 * jl1, k - 2 * kl1, m)) / 12.;
                            df2 = (-q(i + 2 * il1, j + 2 * jl1, k + 2 * kl1, m) + 16. * q(i + il1, j + jl1, k + kl1, m) - 30. * q(i, j, k, m) + 16. * q(i - il1, j - jl1, k - kl1, m) - q(i - 2 * il1, j - 2 * jl1, k - 2 * kl1, m)) / 12.;
                            df3 = (-q(i + 2 * il1, j + 2 * jl1, k + 2 * kl1, m) + 2. * q(i + il1, j + jl1, k + kl1, m) - 2. * q(i - il1, j - jl1, k - kl1, m) + q(i - 2 * il1, j - 2 * jl1, k - 2 * kl1, m)) / 2.;
                            df4 = (q(i + 2 * il1, j + 2 * jl1, k + 2 * kl1, m) - 4. * q(i + il1, j + jl1, k + kl1, m) + 6. * q(i, j, k, m) - 4. * q(i - il1, j - jl1, k - kl1, m) + q(i - 2 * il1, j - 2 * jl1, k - 2 * kl1, m));
                            kn = sqrt((abs(df3) + abs(df4)) / (abs(df1) + abs(df2) + 1.0e-3));
                            psi = 1. - abs(beta20(i, j, k, m) - beta21(i, j, k, m)) / (beta20(i, j, k, m) + beta21(i, j, k, m) + EPSILON);
                            psi = min(1.0, psi / psic);
                            kmkc = kc - kn;
                            signkmkc = 2.0 * (((signed char*)&kmkc)[sizeof(kmkc) - 1] >> 7) + 1.0;
                            sigma = (1.0 + signkmkc) / 2. + (1.0 - signkmkc) / 2. * psi;
                            cbeta0(i, j, k, m) = (cBeta0Min + sigma * sigma * (cBeta0Max - cBeta0Min));

                        }
                    }
                }
            }

            //! Smooth coefficient.
            beta0(I, J, K, M) = beta20(I, J, K, M) + cbeta0 * (q(I + il1, J + jl1, K + kl1, M) - 2.0 * q(I, J, K, M) + q(I - il1, J - jl1, K - kl1, M)) * (q(I + il1, J + jl1, K + kl1, M) - 2.0 * q(I, J, K, M) + q(I - il1, J - jl1, K - kl1, M));
            beta1(I, J, K, M) = beta21(I, J, K, M) + cBeta1 * (q(I, J, K, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M)) * (q(I, J, K, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M));

            //! Global coefficient.
            tau(I, J, K, M) = ((2. * q(I - il1, J - jl1, K - kl1, M) - 3. * q(I, J, K, M) + 1. * q(I + il1, J + jl1, K + kl1, M)) * (q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M) - 3. * q(I - il1, J - jl1, K - kl1, M) + 3. * q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M)));

            //! tau/beta
            w0(I, J, K, M) = tau / (EPSILON + beta0);
            w1(I, J, K, M) = tau / (EPSILON + beta1);

            //! Non-linear weights(not normalized)
            alpha0(I, J, K, M) = DK02 * (1.0 + cAlpha * w0(I, J, K, M) * w0(I, J, K, M));
            alpha1(I, J, K, M) = DK12 * (1.0 + cAlpha * w1(I, J, K, M) * w1(I, J, K, M));

            //! Candidate template.
            q02(I, J, K, M) = -1.0 / 2.0 * q(I + il1, J + jl1, K + kl1, M) + 3.0 / 2.0 * q(I, J, K, M);
            q12(I, J, K, M) = 1.0 / 2.0 * q(I, J, K, M) + 1.0 / 2.0 * q(I - il1, J - jl1, K - kl1, M);
            //! Compute ql.
            qr(I, J, K, M) = (alpha0(I, J, K, M) * q02(I, J, K, M) + alpha1(I, J, K, M) * q12(I, J, K, M)) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));
        }
        else if (structureLimiter == ILMT_WENN3_ZES3)
        {
            //! Define WENN3_ZES3 constants.
            const RDouble EPSILON = 1.0e-40;
            const RDouble DK02 = 1.0 / 3.0, DK12 = 2.0 / 3.0;
            const RDouble cBeta0 = 0.5, cBeta1 = 0.15, cAlpha = 0.4;

            RDouble4D q02(I, J, K, M, fortranArray), q12(I, J, K, M, fortranArray), alpha0(I, J, K, M, fortranArray), alpha1(I, J, K, M, fortranArray);
            RDouble4D tau(I, J, K, M, fortranArray), beta0(I, J, K, M, fortranArray), beta1(I, J, K, M, fortranArray), w0(I, J, K, M, fortranArray), w1(I, J, K, M, fortranArray), tmp(I, J, K, M, fortranArray);

            //! Smooth coefficient.
            beta0(I, J, K, M) = (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)) * (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)) + cBeta0 * (q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I, J, K, M)) * (q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I, J, K, M));
            beta1(I, J, K, M) = (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M)) * (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M)) + cBeta1 * (q(I, J, K, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M)) * (q(I, J, K, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M));

            //! Global smooth coefficient.
            tau(I, J, K, M) = (2. * q(I + il1, J + jl1, K + kl1, M) - 3. * q(I, J, K, M) + 1. * q(I - il1, J - jl1, K - kl1, M)) * (q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M) - 3. * q(I + il1, J + jl1, K + kl1, M) + 3. * q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M));

            //! tau(I, J, K, M)/beta
            w0(I, J, K, M) = tau(I, J, K, M) / (EPSILON + beta0(I, J, K, M));
            w1(I, J, K, M) = tau(I, J, K, M) / (EPSILON + beta1(I, J, K, M));

            //! Non-linear weights(not normalized).
            alpha0(I, J, K, M) = DK02 * (1.0 + cAlpha * w0(I, J, K, M) * w0(I, J, K, M));
            alpha1(I, J, K, M) = DK12 * (1.0 + cAlpha * w1(I, J, K, M) * w1(I, J, K, M));

            //! Candidate template. 
            q02(I, J, K, M) = -1.0 / 2.0 * q(I - il1, J - jl1, K - kl1, M) + 3.0 / 2.0 * q(I, J, K, M);
            q12(I, J, K, M) = 1.0 / 2.0 * q(I, J, K, M) + 1.0 / 2.0 * q(I + il1, J + jl1, K + kl1, M);

            //! Compute ql and normalized.
            ql(I, J, K, M) = (alpha0(I, J, K, M) * q02(I, J, K, M) + alpha1(I, J, K, M) * q12(I, J, K, M)) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));

            //! Smooth coefficient.
            beta0(I, J, K, M) = (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M)) * (q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M)) + cBeta0 * (q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I, J, K, M)) * (q(I + 2 * il1, J + 2 * jl1, K + 2 * kl1, M) - 2.0 * q(I + il1, J + jl1, K + kl1, M) + q(I, J, K, M));
            beta1(I, J, K, M) = (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)) * (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)) + cBeta1 * (q(I, J, K, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M)) * (q(I, J, K, M) - 2.0 * q(I - il1, J - jl1, K - kl1, M) + q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M));

            //! Global smooth coefficient.
            tau(I, J, K, M) = 1. * (2. * q(I - il1, J - jl1, K - kl1, M) - 3. * q(I, J, K, M) + 1. * q(I + il1, J + jl1, K + kl1, M)) * (q(I - 2 * il1, J - 2 * jl1, K - 2 * kl1, M) - 3. * q(I - il1, J - jl1, K - kl1, M) + 3. * q(I, J, K, M) - q(I + il1, J + jl1, K + kl1, M));

            //! tau/beta
            w0(I, J, K, M) = (tau(I, J, K, M) / (EPSILON + beta0(I, J, K, M)));
            w1(I, J, K, M) = (tau(I, J, K, M) / (EPSILON + beta1(I, J, K, M)));

            //! Non-linear weights(not normalized)
            alpha0(I, J, K, M) = DK02 * (1.0 + cAlpha * w0(I, J, K, M) * w0(I, J, K, M));
            alpha1(I, J, K, M) = DK12 * (1.0 + cAlpha * w1(I, J, K, M) * w1(I, J, K, M));

            //! Candidate template.
            q02(I, J, K, M) = -1.0 / 2.0 * q(I + il1, J + jl1, K + kl1, M) + 3.0 / 2.0 * q(I, J, K, M);
            q12(I, J, K, M) = 1.0 / 2.0 * q(I, J, K, M) + 1.0 / 2.0 * q(I - il1, J - jl1, K - kl1, M);

            //! Compute qr and normalized.
            qr(I, J, K, M) = (alpha0(I, J, K, M) * q02(I, J, K, M) + alpha1(I, J, K, M) * q12(I, J, K, M)) / (alpha0(I, J, K, M) + alpha1(I, J, K, M));
        }
        else if (structureLimiter == ILMT_3rd_Minmod_SMOOTH)
        {
            dql(I, J, K, M) = q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M);
            dqr(I, J, K, M) = q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M);

            ql(I, J, K, M) = q(I, J, K, M) + half * Third_Order_Minmod_Smooth(dqr(I, J, K, M), dql(I, J, K, M));
            qr(I, J, K, M) = q(I, J, K, M) - half * Third_Order_Minmod_Smooth(dql(I, J, K, M), dqr(I, J, K, M));
        }
        else if (structureLimiter == ILMT_SMOOTH)
        {
            dql(I, J, K, M) = q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M);
            dqr(I, J, K, M) = q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M);

            dqlr(I, J, K, M) = Smooth(dql(I, J, K, M), dqr(I, J, K, M));

            ql(I, J, K, M) = q(I, J, K, M) + 0.25 * dqlr(I, J, K, M) * ((1.0 - MUSCLCoefXk * dqlr(I, J, K, M)) * dql(I, J, K, M) + (1.0 + MUSCLCoefXk * dqlr(I, J, K, M)) * dqr(I, J, K, M));
            qr(I, J, K, M) = q(I, J, K, M) - 0.25 * dqlr(I, J, K, M) * ((1.0 - MUSCLCoefXk * dqlr(I, J, K, M)) * dqr(I, J, K, M) + (1.0 + MUSCLCoefXk * dqlr(I, J, K, M)) * dql(I, J, K, M));
        }
        else
        {
            if (structureLimiter == ILMT_MINMOD)
            {
                dql(I, J, K, M) = MinMod(q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M), MUSCLCoefXb * (q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M)));
                dqr(I, J, K, M) = MinMod(q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M), MUSCLCoefXb * (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)));
            }
            else if (structureLimiter == ILMT_VAN_ALBADA)
            {
                dql(I, J, K, M) = Vanalbada(q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M), q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M));
                dqr(I, J, K, M) = Vanalbada(q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M), q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M));
            }
            else if (structureLimiter == ILMT_VAN_ALBADA_CLZ)
            {
                dql(I, J, K, M) = VanalbadaCLZ(q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M), q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M));
                dqr(I, J, K, M) = VanalbadaCLZ(q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M), q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M));
            }
            else if (structureLimiter == ILMT_VANLEER)
            {
                dql(I, J, K, M) = Vanleer(q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M), MUSCLCoefXb * (q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M)));
                dqr(I, J, K, M) = Vanleer(q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M), MUSCLCoefXb * (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)));
            }
            else if (structureLimiter == ILMT_MIN_VAN)
            {
                dql(I, J, K, M) = Minvan(q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M), MUSCLCoefXb * (q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M)));
                dqr(I, J, K, M) = Minvan(q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M), MUSCLCoefXb * (q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M)));
            }
            //clz end
            else if (structureLimiter == ILMT_NOLIM)
            {
                dql(I, J, K, M) = q(I, J, K, M) - q(I - il1, J - jl1, K - kl1, M);
                dqr(I, J, K, M) = q(I + il1, J + jl1, K + kl1, M) - q(I, J, K, M);
            }
            else if (structureLimiter == ILMT_FIRST)
            {
                dql(I, J, K, M) = 0;
                dqr(I, J, K, M) = 0;
            }

            ql(I, J, K, M) = q(I, J, K, M) + c1 * dql(I, J, K, M) + c2 * dqr(I, J, K, M);
            qr(I, J, K, M) = q(I, J, K, M) - c2 * dql(I, J, K, M) - c1 * dqr(I, J, K, M);
        }
        delete qProxy;
    }

    void NSSolverStruct::ComputeInviscidFlux(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface)
    {
        StructGrid *grid = StructGridCast(gridIn);

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nInviscidFluxModify = parameters->GetnInviscidFluxModify();

        FieldProxy *fluxProxy = CreateFieldProxy(grid);

        int structureScheme = GlobalDataBase::GetIntParaFromDB("str_scheme");
        INVScheme inviscidScheme = GetINVScheme(structureScheme);

        //! To compute the inviscid flux of the faces using the QL and QR.
        InviscidFluxWrap(grid, qlProxy, qrProxy, fluxProxy, inviscidScheme, iSurface);

#ifdef USE_HYFLOW_CODE
        if (nInviscidFluxModify == 1)
            ModifyFlux(gridIn, fluxProxy, iSurface);   //by dms
#endif
//! To compute the summation of fluxes on faces.
        LoadFlux(grid, fluxProxy, iSurface);

        delete fluxProxy;
    }

    void NSSolverStruct::CorrectInviscidFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
    {
        Param_NSSolverStruct *parameters = GetControlParameters();
        int nChemical = parameters->GetChemicalFlag();
        if (nChemical == 0)
        {//! Apply the original method, it is no need to modify the flux on wall.
            return;
        }
        RDouble wallTemperature = parameters->GetWallTemperature();
        if (wallTemperature < 0.0)
        {//! The adiabatic wall applies the original method.
            return;
        }

        int viscousType = parameters->GetViscousType();
        RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
        //! The non-dimensional value of wall temperature.
        RDouble nondimensionalTemperatureWall = wallTemperature / refDimensionalTemperature;
        RDouble temperatureWall[3] = { nondimensionalTemperatureWall, nondimensionalTemperatureWall, nondimensionalTemperatureWall };

        //! Otherwise, to modify the flux of solid boundaries.
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
        RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
        RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());
        RDouble4D &faceArea = *(grid->GetFaceArea());
        RDouble4D &faceNormalVelocity = *(grid->GetFaceNormalVelocity());

        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &faceFluxes = fluxProxy->GetField_STR();

        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
        int nSpecies = nLaminar + nChemical - nNSEquation;

        RDouble catalyticCoef = parameters->GetCatalyticCoef();
        RDouble3D *surfaceMassFraction = 0;
        if (nChemical > 0 && nSpecies > 0)    //! The fully catalytic wall condition.
        {
            surfaceMassFraction = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("surfaceMassFraction"));
        }
        int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
        RDouble3D *surfaceSlipVariables = NULL;
        if (nSlipBCModel > 0)
        {
            surfaceSlipVariables = reinterpret_cast <RDouble3D *> (gridIn->GetDataPtr("surfaceSlipVariables"));
        }

        //! Save the variables and inviscid fluxes on wall.
        RDouble *wallVariables = new RDouble[nEquation]();
        RDouble *wallFluxes = new RDouble[nEquation]();

        using namespace IDX;
        RDouble gasConstant = gas->GetUniversalGasConstant();

        RDouble normalComponentX = 0.0, normalComponentY = 0.0, normalComponentZ = 0.0, surfaceArea = 0.0;
        RDouble relativeVelocity = 0.0, faceMotionVelocity = 0.0, totalEnthalpy = 0.0, vibrationEnergy = 0.0, electronEnergy = 0.0;

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        int indexOfWall = 0, indexOfCell = 0;

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int nBCType = structBC->GetBCType();

            if (!(IsWall(nBCType) && viscousType > INVISCID) && nBCType != PHENGLEI::ABLATION_SURFACE)
            {
                continue;
            }

            int *faceDirectionIndex = structBC->GetFaceDirectionIndex();
            //! 0 stands for I-direction, 1 stands for J-direction, and 2 for K-direction.
            int directionIndex = structBC->GetFaceDirection() + 1;
            //! Judge the face is on the left side or right side in the current zone,-1 stands for left side, and 1 for right side.

            int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            if (directionIndex != iSurface)
            {
                continue;
            }

            int id = 0, jd = 0, kd = 0;
            GetBCFaceIDX(faceDirectionIndex, id, jd, kd);

            kst = kst + kd;
            ked = ked + kd;
            jst = jst + jd;
            jed = jed + jd;
            ist = ist + id;
            ied = ied + id;

            indexOfCell = 0;
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        //! To obtain the variables on wall.
                        if (nSlipBCModel > 0)
                        {
                            for (int m = 0; m < 3; ++ m)
                            {
                                temperatureWall[m] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m);
                                wallVariables[m + 1] = (*surfaceSlipVariables)(indexOfWall, indexOfCell, m + 3);
                            }
                        }
                        else
                        {
                            if (viscousType > INVISCID)
                            {
                                wallVariables[IU] = 0.0;
                                wallVariables[IV] = 0.0;
                                wallVariables[IW] = 0.0;
                            }
                            else
                            {
                                wallVariables[IU] = primitiveVariables(i, j, k, IU);
                                wallVariables[IV] = primitiveVariables(i, j, k, IV);
                                wallVariables[IW] = primitiveVariables(i, j, k, IW);
                            }
                        }
                        wallVariables[IP] = primitiveVariables(i, j, k, IP);

                        if (ABS(catalyticCoef) <= EPSILON && nBCType != PHENGLEI::ABLATION_SURFACE) //Fully non-catalytic wall.
                        {
                            for (int m = nNSEquation; m < nNSEquation + nSpecies; ++ m)
                            {
                                wallVariables[m] = primitiveVariables(i, j, k, m);
                            }
                        }
                        else //Fully catalytic wall or finite catalytic wall.
                        {
                            for (int m = nNSEquation; m < nNSEquation + nSpecies; ++ m)
                            {
                                wallVariables[m] = (*surfaceMassFraction)(indexOfWall, indexOfCell, m - nNSEquation);
                            }
                        }

                        //! The isothermal condition reveals that both vibration temperature and electron temperature \n.
                        //! equal to translation-rotation temperature. Tve=Tv=Te=Tw.

                        //! Compute density on wall.
                        RDouble massReciprocal = gas->ComputeMolecularWeightReciprocal(wallVariables);
                        wallVariables[IR] = wallVariables[IP] / (gasConstant * temperatureWall[ITT] * massReciprocal);

                        //! The geometric variables on wall.
                        normalComponentX = faceNormalComponentX(i, j, k, iSurface);
                        normalComponentY = faceNormalComponentY(i, j, k, iSurface);
                        normalComponentZ = faceNormalComponentZ(i, j, k, iSurface);
                        surfaceArea = faceArea(i, j, k, iSurface);
                        faceMotionVelocity = faceNormalVelocity(i, j, k, iSurface);

                        //! Obtain the velocity on wall.
                        relativeVelocity = wallVariables[IU] * normalComponentX + wallVariables[IV] * normalComponentY + wallVariables[IW] * normalComponentZ - faceMotionVelocity;
                        //! Compute the total enthalpy.
                        totalEnthalpy = gas->GetMixedGasEnthalpy(wallVariables, temperatureWall[ITT], temperatureWall[ITV], temperatureWall[ITE]); //Tve=Tv=Te=Tw
                        totalEnthalpy += 0.5 * (pow(wallVariables[IU], 2) + pow(wallVariables[IV], 2) + pow(wallVariables[IW], 2));

                        wallFluxes[IR] = wallVariables[IR] * relativeVelocity;
                        wallFluxes[IU] = wallFluxes[IR] * wallVariables[IU] + wallVariables[IP] * normalComponentX;
                        wallFluxes[IV] = wallFluxes[IR] * wallVariables[IV] + wallVariables[IP] * normalComponentY;
                        wallFluxes[IW] = wallFluxes[IR] * wallVariables[IW] + wallVariables[IP] * normalComponentZ;
                        wallFluxes[IP] = wallFluxes[IR] * totalEnthalpy + wallVariables[IP] * faceMotionVelocity;
                        for (int m = nNSEquation; m < nLaminar; ++ m)
                        {
                            wallFluxes[m] = wallFluxes[IR] * wallVariables[m];
                        }

                        if (nChemical == 1)
                        {
                            wallFluxes[nLaminar] = 0.0;
                        }
                        if (nTemperatureModel > 1)
                        {
                            //! Compute the internal energy flux of vibration and electron.
                            vibrationEnergy = gas->GetMixedGasVibrationEnergy(wallVariables, temperatureWall[ITV]);   //! Tve=Tv=Te=Tw.
                            electronEnergy = gas->GetMixedGasElectronEnergy(wallVariables, temperatureWall[ITE]);     //! Tve=Tv=Te=Tw.

                            RDouble ceDivideMe = 0.0, peDivideDensity = 0.0;
                            massReciprocal = gas->ComputeMolecularWeightReciprocalWithoutElectron(wallVariables, ceDivideMe);
                            peDivideDensity = gasConstant * temperatureWall[ITE] * ceDivideMe;     //! Tve=Tv=Te=Tw.

                            if (nTemperatureModel == 2)         //! Two-Temperature model.
                            {
                                wallFluxes[nLaminar + nChemical] = wallFluxes[IR] * (vibrationEnergy + electronEnergy + peDivideDensity) + wallVariables[IR] * peDivideDensity * faceMotionVelocity;
                            }
                            else if (nTemperatureModel == 3)    //! Three-Temperature model.
                            {
                                wallFluxes[nLaminar + nChemical] = wallFluxes[IR] * vibrationEnergy;
                                wallFluxes[nLaminar + nChemical + 1] = wallFluxes[IR] * (electronEnergy + peDivideDensity) + wallVariables[IR] * peDivideDensity * faceMotionVelocity;
                            }
                        }

                        //! Modify inviscid flux.
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            faceFluxes(i, j, k, m) = wallFluxes[m] * surfaceArea;
#ifdef USE_ErrorDataDebug
                            if (isErrorData(faceFluxes(i, j, k, m)))
                            {
                                int zoneID = grid->GetZoneID();
                                PrintToWindow(" zoneID = ", zoneID, "i = ", i, "j = ", j, "k = ", k, "m = ", m);
                                TK_Exit::ExceptionExit("This is in function CorrectInviscidFlux", true);
                            }
#endif
                        }
                        //Next grid cell.
                        ++ indexOfCell;
                    }//i end.
                } //j end.
            } //k end.
            ++ indexOfWall; //Next surface regions.
        }//iBCRegion end.
        delete[]wallVariables;    wallVariables = NULL;
        delete[]wallFluxes;    wallFluxes = NULL;
    }

    FaceProxy *NSSolverStruct::CreateFaceProxy(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int nEquation = GetNumberOfEquations();

        int nlen = MAX(MAX(ni, nj), nk);

        FaceProxy *faceProxy = new FaceProxy();
        faceProxy->Create(nlen, nEquation);

        return faceProxy;
    }

    GeomProxy *NSSolverStruct::CreateGeomProxy(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int nlen = MAX(MAX(ni, nj), nk);

        GeomProxy *geomProxy = new GeomProxy();
        geomProxy->Create(nlen);

        return geomProxy;
    }

    bool iterijk(int &i, int &j, int &k, int ist, int ied, int jst, int jed, int kst, int ked)
    {
        ++ i;
        if (i <= ied)
        {
            return true;
        }
        else
        {
            i = ist;
            ++ j;
            if (j <= jed)
            {
                return true;
            }
            else
            {
                j = jst;
                ++ k;
                if (k <= ked)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
        }
    }

    bool NSSolverStruct::GetInviscidFaceValue(Grid *gridIn, FaceProxy *faceProxy, StructIter *structIter, FieldProxy *qlProxy, FieldProxy *qrProxy, int nlen)
    {
        Param_NSSolverStruct *parameters = GetControlParameters();
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &ql = qlProxy->GetField_STR();
        RDouble4D &qr = qrProxy->GetField_STR();

        RDouble4D &gxfn = *(grid->GetFaceNormalX());
        RDouble4D &gyfn = *(grid->GetFaceNormalY());
        RDouble4D &gzfn = *(grid->GetFaceNormalZ());
        RDouble4D &garea = *(grid->GetFaceArea());
        RDouble4D &gvgn = *(grid->GetFaceNormalVelocity());

        int nEquation = GetNumberOfEquations();
        int nEnergyRecycle = parameters->GetnEnergyRecycle();
        int ntmodel = parameters->GetTemperatureModel();
        int nLaminar = parameters->GetLaminarNumber();
        RDouble densityMin = parameters->GetDensityMin();

        int nQlLimitMethod = parameters->GetnQlLimitMethod();

        RDouble3D &gamma = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));
        RDouble3D &rtem = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("rtem"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

        int mTT = parameters->GetmTT();
        int mTV = parameters->GetmTV();
        int mTE = parameters->GetmTE();
        using namespace IDX;
        int nm = parameters->GetNSEquationNumber();

        GeomProxy *geomProxy = faceProxy->GetGeomProxy();

        RDouble *xfn = geomProxy->GetFaceNormalX();
        RDouble *yfn = geomProxy->GetFaceNormalY();
        RDouble *zfn = geomProxy->GetFaceNormalZ();
        RDouble *area = geomProxy->GetFaceArea();
        RDouble *vgn = geomProxy->GetFaceVelocity();

        RDouble *gamal = faceProxy->GetGamaL();
        RDouble *gamar = faceProxy->GetGamaR();

        RDouble **temperatureL = faceProxy->GetLeftTemperature();
        RDouble **temperatureR = faceProxy->GetRightTemperature();

        RDouble **face_ql = faceProxy->GetQL();
        RDouble **face_qr = faceProxy->GetQR();

        RDouble *lrtl = faceProxy->GetPressureCoefficientL();
        RDouble *lrtr = faceProxy->GetPressureCoefficientR();
        RDouble *primitiveVarL = new RDouble[nEquation];
        RDouble *primitiveVarR = new RDouble[nEquation];
        RDouble leftTemperature[3], rightTemperature[3];

        bool flag = false;

        int i = structIter->i0;
        int j = structIter->j0;
        int k = structIter->k0;

        int il1 = structIter->il1;
        int jl1 = structIter->jl1;
        int kl1 = structIter->kl1;

        int iSurface = structIter->nsurf;

        int ist = structIter->ist;
        int ied = structIter->ied;
        int jst = structIter->jst;
        int jed = structIter->jed;
        int kst = structIter->kst;
        int ked = structIter->ked;

        int ix = 0;

        do
        {
            int il = i + il1;
            int jl = j + jl1;
            int kl = k + kl1;

            #ifdef USE_HYFLOW_CODE
            if (nEnergyRecycle > 0)
            {
                if (nQlLimitMethod == 0)
                {
                    if (ql(i, j, k, IR) <= 0.0) ql(i, j, k, IR) = q(i, j, k, IR);
                    if (ql(i, j, k, IP) <= 0.0) ql(i, j, k, IP) = q(i, j, k, IP);

                    if (qr(il, jl, kl, IR) <= 0.0) qr(il, jl, kl, IR) = q(il, jl, kl, IR);
                    if (qr(il, jl, kl, IP) <= 0.0) qr(il, jl, kl, IP) = q(il, jl, kl, IP);

                    for (int n = nm; n < nEquation; ++ n)
                    {
                        if (ql(i, j, k, n) <= 0.0) ql(i, j, k, n) = q(i, j, k, n);
                        if (qr(il, jl, kl, n) <= 0.0) qr(il, jl, kl, n) = q(il, jl, kl, n);
                    }
                }
                else if (nQlLimitMethod == 1)
                {
                    if (ql(i, j, k, IR) <= 0.0) ql(i, j, k, IR) = densityMin;
                    if (ql(i, j, k, IP) <= 0.0) ql(i, j, k, IP) = 1.0e-25;

                    if (qr(il, jl, kl, IR) <= 0.0) qr(il, jl, kl, IR) = densityMin;
                    if (qr(il, jl, kl, IP) <= 0.0) qr(il, jl, kl, IP) = 1.0e-25;

                    for (int n = nm; n < nEquation; ++ n)
                    {
                        if (ql(i, j, k, n) <= 0.0) ql(i, j, k, n) = 1.0e-25;
                        if (qr(il, jl, kl, n) <= 0.0) qr(il, jl, kl, n) = 1.0e-25;
                    }
                }
                else if (nQlLimitMethod == 2)
                {
                    ModifyFluxRL(ql(i, j, k, IR), qr(il, jl, kl, IR), q(i, j, k, IR), q(il, jl, kl, IR));
                    ModifyFluxRL(ql(i, j, k, IP), qr(il, jl, kl, IP), q(i, j, k, IP), q(il, jl, kl, IP));

                    for (int n = nm; n < nEquation; ++ n)
                    {
                        ModifyFluxRL(ql(i, j, k, n), qr(il, jl, kl, n), q(i, j, k, n), q(il, jl, kl, n));
                    }
                }
            }
            #endif

            for (int m = 0; m < nEquation; ++ m)
            {
                face_ql[m][ix] = ql(i, j, k, m);
                face_qr[m][ix] = qr(il, jl, kl, m);

                primitiveVarL[m] = face_ql[m][ix];
                primitiveVarR[m] = face_qr[m][ix];
            }

            lrtl[ix] = rtem(i, j, k);
            lrtr[ix] = rtem(il, jl, kl);

            gamal[ix] = gamma(i, j, k);
            gamar[ix] = gamma(il, jl, kl);

            leftTemperature[0] = t(i, j, k, mTT);
            rightTemperature[0] = t(il, jl, kl, mTT);

            leftTemperature[1] = t(i, j, k, mTV);
            rightTemperature[1] = t(il, jl, kl, mTV);

            leftTemperature[2] = t(i, j, k, mTE);
            rightTemperature[2] = t(il, jl, kl, mTE);

            //! To obtain the left and right temperature.
            if (nEnergyRecycle == 1)
            {
                gas->GetTemperatureR(primitiveVarL, leftTemperature[0], leftTemperature[1], leftTemperature[2]);
                gas->GetTemperatureR(primitiveVarR, rightTemperature[0], rightTemperature[1], rightTemperature[2]);
            }
            else if (nEnergyRecycle == 0)
            {
                gas->GetTemperature(primitiveVarL, leftTemperature[0], leftTemperature[1], leftTemperature[2]);
                gas->GetTemperature(primitiveVarR, rightTemperature[0], rightTemperature[1], rightTemperature[2]);
            }
            /*else
            {
                gas->GetTemperatureR2(primitiveVarL, leftTemperature[0], leftTemperature[1], leftTemperature[2]);
                gas->GetTemperatureR2(primitiveVarR, rightTemperature[0], rightTemperature[1], rightTemperature[2]);

                for (int m = 0; m < nEquation; ++  m)
                {
                    face_ql[m][ix] = primitiveVarL[m];
                    face_qr[m][ix] = primitiveVarR[m];
                }
            }*/


            xfn[ix] = gxfn(il, jl, kl, iSurface);
            yfn[ix] = gyfn(il, jl, kl, iSurface);
            zfn[ix] = gzfn(il, jl, kl, iSurface);
            area[ix] = garea(il, jl, kl, iSurface);
            vgn[ix] = gvgn(il, jl, kl, iSurface);

            temperatureL[0][ix] = leftTemperature[mTT];
            temperatureR[0][ix] = rightTemperature[mTT];

            temperatureL[1][ix] = leftTemperature[mTV];
            temperatureR[1][ix] = rightTemperature[mTV];

            temperatureL[2][ix] = leftTemperature[mTE];
            temperatureR[2][ix] = rightTemperature[mTE];

            ++ ix;

            flag = iterijk(i, j, k, ist, ied, jst, jed, kst, ked);
        } while (flag && ix < nlen);

        delete[] primitiveVarL;
        delete[] primitiveVarR;

        structIter->i1 = i;
        structIter->j1 = j;
        structIter->k1 = k;

        structIter->nsize = ix;

        return flag;
    }

    bool NSSolverStruct::LoadFaceFlux(Grid *gridIn, FaceProxy *faceProxy, StructIter *structIter, FieldProxy *fluxProxy, int nlen)
    {
        RDouble4D &flux = fluxProxy->GetField_STR();

        int nEquation = GetNumberOfEquations();

        GeomProxy *geomProxy = faceProxy->GetGeomProxy();
        RDouble *area = geomProxy->GetFaceArea();
        RDouble **faceFlux = faceProxy->GetFlux();

        bool flag = false;

        int i = structIter->i0;
        int j = structIter->j0;
        int k = structIter->k0;

        int il1 = structIter->il1;
        int jl1 = structIter->jl1;
        int kl1 = structIter->kl1;

        int ist = structIter->ist;
        int ied = structIter->ied;
        int jst = structIter->jst;
        int jed = structIter->jed;
        int kst = structIter->kst;
        int ked = structIter->ked;

        int il, jl, kl;

        int ix = 0;

        do
        {
            il = i + il1;
            jl = j + jl1;
            kl = k + kl1;

            RDouble &areas = area[ix];

            for (int m = 0; m < nEquation; ++ m)
            {
                flux(il, jl, kl, m) = faceFlux[m][ix] * areas;
            }
            ++ ix;

            flag = iterijk(i, j, k, ist, ied, jst, jed, kst, ked);
        } while (flag && ix < nlen);

        structIter->i1 = i;
        structIter->j1 = j;
        structIter->k1 = k;

        structIter->nsize = ix;

        return flag;
    }

    void NSSolverStruct::InviscidFluxWrap(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, FieldProxy *fluxProxy, INVScheme invScheme, int iSurface)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        FaceProxy *faceProxy = CreateFaceProxy(grid);
        faceProxy->SetGeomProxy(CreateGeomProxy(grid));

        Param_NSSolverStruct *parameters = this->GetControlParameters();

        int nlen = MAX(MAX(ni, nj), nk);
        //int nlen = 10;

        StructIter *structIter = CreateInvStructIter(grid, iSurface);

        int nface = 0;

        do
        {
            bool flag1 = GetInviscidFaceValue(grid, faceProxy, structIter, qlProxy, qrProxy, nlen);

            faceProxy->setsize(structIter->nsize);
            nface += structIter->nsize;

            InviscidSchemeParameter invSchemeParaProxy;
            SetInviscidSchemeParameters(&invSchemeParaProxy, faceProxy, parameters);

            //! Compute the fluxes on the faces using the QL and QR.
            //! The variable faceProxy stores the new computation values.
            invScheme(faceProxy, &invSchemeParaProxy);

            bool hasNext = LoadFaceFlux(grid, faceProxy, structIter, fluxProxy, nlen);

            if (!hasNext) break;

            structIter->i0 = structIter->i1;
            structIter->j0 = structIter->j1;
            structIter->k0 = structIter->k1;
        } while (true);

        delete faceProxy;
        delete structIter;
    }

    StructIter *NSSolverStruct::CreateInvStructIter(Grid *gridIn, int iSurface)
    {
        StructGrid *grid = StructGridCast(gridIn);

        int il1, jl1, kl1;
        grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

        int iFaceStart, iFaceEnd, jFaceStart, jFaceEnd, kFaceStart, kFaceEnd;
        grid->GetFaceIterationIndex(iFaceStart, iFaceEnd, jFaceStart, jFaceEnd, kFaceStart, kFaceEnd, iSurface);

        StructIter *structIter = new StructIter();

        structIter->i0 = iFaceStart;
        structIter->j0 = jFaceStart;
        structIter->k0 = kFaceStart;

        structIter->il1 = il1;
        structIter->jl1 = jl1;
        structIter->kl1 = kl1;

        structIter->ist = iFaceStart;
        structIter->ied = iFaceEnd;
        structIter->jst = jFaceStart;
        structIter->jed = jFaceEnd;
        structIter->kst = kFaceStart;
        structIter->ked = kFaceEnd;

        structIter->nsurf = iSurface;

        return structIter;
    }

    void NSSolverStruct::Diagonal(Grid *gridIn)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        RDouble3D &vol = *(grid->GetCellVolume());

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
        int nChemicalRadius = parameters->GetNChemicalRadius();
        int nChemicalSource = parameters->GetNChemicalSource();
        int numberOfSpecies = parameters->GetNumberOfSpecies();

        RDouble ChemicalSpectrumRadiusCoef = parameters->GetChemicalSpectrumRadiusCoef();
        RDouble ViscousSpectrumRadiusCoef = parameters->GetViscousSpectrumRadiusCoef();
        RDouble InviscidSpectrumRadiusCoef = parameters->GetInviscidSpectrumRadiusCoef();

        RDouble4D &diagonal = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
        RDouble4D &diagonal0 = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal0"));
        RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
        RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));

        RDouble4D &srs = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("srs"));
        RDouble3D &dt = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));

        Range M(0, nEquation - 1);

        int nDim = GetDim();

        RDouble odt, rad, rad_euler, rad_ns;
        RDouble wmig = 1.0;
        RDouble beta = 1.0;
        RDouble betw = beta * wmig;

        RDouble dualTimeSpectrumC1 = zero;
        RDouble dualTimeSpectrumC2 = one;

        int isUnsteady = parameters->GetIsUnsteady();
        if (isUnsteady)
        {
            RDouble dualTimeCoefficient[7];
            const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
            ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
            dualTimeSpectrumC1 = -dualTimeCoefficient[3];
            dualTimeSpectrumC2 = dualTimeCoefficient[6];
        }

        RDouble viscousSpectrumC1 = zero;

        int viscousType = parameters->GetViscousType();
        if (viscousType > INVISCID)
        {
            viscousSpectrumC1 = ViscousSpectrumRadiusCoef;
            //      viscousSpectrumC1 = 1.0;                                                                  ///
        }

        //   InviscidSpectrumRadiusCoef = max(0.0, InviscidSpectrumRadiusCoef - 1.0);                       ///

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    odt = dualTimeSpectrumC1 * vol(i, j, k) + dualTimeSpectrumC2 / dt(i, j, k);

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        diagonal(i, j, k, m) = odt;
                    }

                    rad_euler = invSpectralRadius(i, j, k, 0);
                    rad_ns = visSpectralRadius(i, j, k, 0);
                    for (int iSurface = 1; iSurface < nDim; ++ iSurface)
                    {
                        rad_euler += invSpectralRadius(i, j, k, iSurface);
                        rad_ns += visSpectralRadius(i, j, k, iSurface);   //! Modified by clz : 2012-7-19.
                    }
                    //     rad = rad_euler + viscousSpectrumC1 * rad_ns;
                    rad = rad_euler * InviscidSpectrumRadiusCoef + viscousSpectrumC1 * rad_ns;
                    rad *= betw;

                    diagonal0(i, j, k, 0) = rad * InviscidSpectrumRadiusCoef;

                    for (int m = 0; m < nEquation; ++ m)       //! Here is nl, not nEquation.
                    {
                        diagonal(i, j, k, m) += rad;
                    }
                }
            }
        }

        //clz begin
        if (nChemical == 1 && nChemicalSource == 1 && nChemicalRadius == 1)
        {
            int nSpeciesEquation = gas->GetSpeciesEquationNumber();
            for (int k = kCellStart; k <= kCellEnd; ++ k)
            {
                for (int j = jCellStart; j <= jCellEnd; ++ j)
                {
                    for (int i = iCellStart; i <= iCellEnd; ++ i)
                    {

                        for (int s = 0; s < nSpeciesEquation; ++ s)
                        {
                            diagonal(i, j, k, nNSEquation + s) -= srs(i, j, k, s) * ChemicalSpectrumRadiusCoef;
                        }
                        for (int it = 0; it < nTemperatureModel - 1; ++ it)
                        {
                            diagonal(i, j, k, nNSEquation + numberOfSpecies + it) -= srs(i, j, k, numberOfSpecies + it) * ChemicalSpectrumRadiusCoef;
                        }

                        //rad = diagonal0(i, j, k, 0);
                        //for (int m = 0; m < nEquation; ++  m)       //! Here is nl, not nEquation.
                        //{
                        //   diagonal(i, j, k, m) += rad;
                        //}
                    }
                }
            }
        }
        //clz end
        GhostCell3D(diagonal, ni, nj, nk, nEquation);
        //grid->GhostCell3DExceptInterface(diagonal, nEquation);

#ifdef USE_HYFLOW_CODE
 //   GAS_SPACE::modifyDiagonal(grid, nEquation);
#endif

    }

#ifdef SMARTARRAY
    //! With faster, smart array optimization.
    void NSSolverStruct::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
    {
        PrintToWindow("Warning: need to check the optimized code, after consider multi-sweep!");
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &dq = dqProxy->GetField_STR();
        RDouble4D &deltaFlux = LUplusDQ->GetField_STR();

        RDouble4D &gxfn = *(grid->GetFaceNormalX());
        RDouble4D &gyfn = *(grid->GetFaceNormalY());
        RDouble4D &gzfn = *(grid->GetFaceNormalZ());
        RDouble4D &gvgn = *(grid->GetFaceNormalVelocity());
        RDouble4D &garea = *(grid->GetFaceArea());
        RDouble3D &gvol = *(grid->GetCellVolume());

        int *iBlank = grid->GetCellTypeContainer();

        Param_NSSolverStruct *parameters = GetControlParameters();

        int nm = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

        RDouble refGama = parameters->GetRefGama();

        int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

        RDouble *nxyz = new RDouble[4];
        RDouble *d_q = new RDouble[nEquation];
        RDouble *prim = new RDouble[nEquation];
        RDouble *q_loc = new RDouble[nEquation];
        RDouble *rhs0 = new RDouble[nEquation];
        RDouble *de = new RDouble[nEquation];
        RDouble *df = new RDouble[nEquation];
        RDouble *dg = new RDouble[nEquation];
        RDouble *dq_pre = new RDouble[nEquation];
        RDouble **MG = NewPointer2 <RDouble>(nEquation, nEquation);
        RDouble *temperature = new RDouble[nTemperatureModel];
        RDouble *q_pre = new RDouble[nEquation];

        RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
        RDouble4D &diagonal = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
        RDouble4D &res = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));

        //clz begin
        RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
        //clz end

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &dt = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));

        RDouble3D &gamma = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));

        Range II, JJ, KK;
        GetRange(ni, nj, nk, -2, 1, II, JJ, KK);
        Range MM(0, nEquation - 1);

        RDouble wmig = 1.0;
        RDouble beta = 1.0;
        RDouble coevis = 0.0;

        Range I, J, K;
        GetRange(ni, nj, nk, 0, -1, I, J, K);

        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        int viscousType = parameters->GetViscousType();
        //! Obtain the viscosity coefficient, Modified by LiPeng on Mar. 15, 2019.
        RDouble3D viscousLaminar, viscousTurbulence;
        if (viscousType > INVISCID)
        {
            viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
            viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
        }
        //! End of the modification.

        int nDim = GetDim();
        int il, jl, kl;

        RDouble *dq_i = new RDouble[nEquation];
        RDouble *dq_j = new RDouble[nEquation];
        RDouble *dq_k = new RDouble[nEquation];

        RDouble *prim_i = new RDouble[nEquation];
        RDouble *prim_j = new RDouble[nEquation];
        RDouble *prim_k = new RDouble[nEquation];

        RDouble *q_loc_i = new RDouble[nEquation];
        RDouble *q_loc_j = new RDouble[nEquation];
        RDouble *q_loc_k = new RDouble[nEquation];

        RDouble *gykb_i = new RDouble[5];
        RDouble *gykb_j = new RDouble[5];
        RDouble *gykb_k = new RDouble[5];

        RDouble **MG_i = NewPointer2 <RDouble>(nEquation, nEquation);
        RDouble **MG_j = NewPointer2 <RDouble>(nEquation, nEquation);
        RDouble **MG_k = NewPointer2 <RDouble>(nEquation, nEquation);


        for (int k = kst; k <= ked; ++ k)
        {
            int km1 = k - 1;
            if (nDim != THREE_D)
            {
                km1 = 1;
            }
            for (int j = jst; j <= jed; ++ j)
            {
                int jm1 = j - 1;
                for (int i = ist; i <= ied; ++ i)
                {
                    int im1 = i - 1;
                    int cellLabel = (ni - 1) * (nj - 1) * (k - 1) + (ni - 1) * (j - 1) + (i - 1);
                    if (iBlank[cellLabel] == IN_HOLE_COLOR)
                    {
                        for (int m = 0; m < nl; ++ m)
                        {
                            dq(i, j, k, m) = 0.0;
                            //rhs0[m] = 0.0;
                        }
                        continue;
                    }
                    for (int m = 0; m < nl; ++ m)
                    {
                        rhs0[m] = 0.0;
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prim_i[m] = q(im1, j, k, m);
                        prim_j[m] = q(i, jm1, k, m);
                    }
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        q_loc_i[m] = q(im1, j, k, m);
                        q_loc_j[m] = q(i, jm1, k, m);
                    }
                    for (int m = 0; m < nl; ++ m)
                    {
                        dq_i[m] = dq(im1, j, k, m);
                        dq_j[m] = dq(i, jm1, k, m);
                    }
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        q_pre[m] = q(i, j, k, m);
                    }
                    for (int m = 0; m < nl; ++ m)
                    {
                        dq_pre[m] = dq(i, j, k, m);
                    }
                    gykb_i[0] = half * (gxfn(i, j, k, 1) + gxfn(im1, j, k, 1));
                    gykb_i[1] = half * (gyfn(i, j, k, 1) + gyfn(im1, j, k, 1));
                    gykb_i[2] = 0.0;
                    gykb_i[3] = half * (garea(i, j, k, 1) + garea(im1, j, k, 1));
                    gykb_i[4] = half * (gvgn(i, j, k, 1) + gvgn(im1, j, k, 1));

                    gykb_j[0] = half * (gxfn(i, j, k, 2) + gxfn(i, jm1, k, 2));
                    gykb_j[1] = half * (gyfn(i, j, k, 2) + gyfn(i, jm1, k, 2));
                    gykb_j[2] = 0.0;
                    gykb_j[3] = half * (garea(i, j, k, 2) + garea(i, jm1, k, 2));
                    gykb_j[4] = half * (gvgn(i, j, k, 2) + gvgn(i, jm1, k, 2));

                    if (nDim == THREE_D)
                    {
                        gykb_i[2] = half * (gzfn(i, j, k, 1) + gzfn(im1, j, k, 1));
                        gykb_j[2] = half * (gzfn(i, j, k, 2) + gzfn(i, jm1, k, 2));
                    }

                    RDouble &ra = invSpectralRadius(im1, j, k, 0);
                    RDouble &rb = invSpectralRadius(i, jm1, k, 1);

                    RDouble gama1 = gamma(im1, j, k);
                    RDouble gama2 = gamma(i, jm1, k);

                    RDouble tm_i = t(im1, j, k, 0);
                    RDouble tm_j = t(i, jm1, k, 0);

                    if (ifLowSpeedPrecon == 0)
                    {
                        //! The previous method of computing the (M + R) * dQ.
                        MXDQ_STD(prim_i, gykb_i, dq_i, de, gama1, nm, nl, ra, 1);
                        MXDQ_STD(prim_j, gykb_j, dq_j, df, gama2, nm, nl, rb, 1);
                    }
                    else
                    {
                        PreconMatrixFreeDF(q_loc_i, dq_i, gykb_i, MG_i, tm_i, refGama, de, nl, 1);
                        PreconMatrixFreeDF(q_loc_j, dq_j, gykb_j, MG_j, tm_j, refGama, df, nl, 1);
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        dg[m] = 0.0;
                    }

                    if (nDim == THREE_D)
                    {
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            prim_k[m] = q(i, j, km1, m);
                        }
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            q_loc_k[m] = q(i, j, km1, m);
                        }
                        for (int m = 0; m < nl; ++ m)
                        {
                            dq_k[m] = dq(i, j, km1, m);
                        }
                        gykb_k[0] = half * (gxfn(i, j, k, 3) + gxfn(i, j, km1, 3));
                        gykb_k[1] = half * (gyfn(i, j, k, 3) + gyfn(i, j, km1, 3));
                        gykb_k[2] = half * (gzfn(i, j, k, 3) + gzfn(i, j, km1, 3));
                        gykb_k[3] = half * (garea(i, j, k, 3) + garea(i, j, km1, 3));
                        gykb_k[4] = half * (gvgn(i, j, k, 3) + gvgn(i, j, km1, 3));
                        RDouble &rc = invSpectralRadius(i, j, km1, 2);
                        RDouble gama3 = gamma(i, j, km1);
                        RDouble tm_k = t(i, j, km1, 0);
                        if (ifLowSpeedPrecon == 0)
                        {
                            MXDQ_STD(prim_k, gykb_k, dq_k, dg, gama3, nm, nl, rc, 1);
                        }
                        else
                        {
                            PreconMatrixFreeDF(q_loc_k, dq_k, gykb_k, MG_k, tm_k, refGama, dg, nl, 1);
                        }
                    }

                    //! Solve rhs0.
                    for (int m = 0; m < nl; ++ m)
                    {
                        rhs0[m] += -wmig * (de[m] + df[m] + dg[m]);  //! The first method.
                    }

                    if (viscousType > INVISCID)         //! Implicit treatment of viscosity term,added by clz : 2012-7-19，modified by LiPeng on Jan 21, 2019.
                    {

                        RDouble rad_vis_i = visSpectralRadius(im1, j, k, 0);    //! The viscous spectrum radius.
                        RDouble rad_vis_j = visSpectralRadius(i, jm1, k, 1);    //! The viscous spectrum radius.
                        RDouble rad_vis_k = 0.0;
                        if (nDim == THREE_D)
                        {
                            rad_vis_k = visSpectralRadius(i, j, km1, 2);
                        }

                        for (int m = 0; m < nl; ++ m)
                        {
                            rhs0[m] += -0.5 * (rad_vis_i * dq_i[m] + rad_vis_j * dq_j[m] + rad_vis_k * dq_k[m]); //the second method
                        }
                    }

                    //! Solve dqa.
                    if (ifLowSpeedPrecon == 1)
                    {
                        RDouble tm = t(i, j, k, 0);
                        ComputePreconMatrix(1, MG, q_pre, tm, refGama);
                        PreconMatrixDF(dq_pre, MG, nl);
                    }

                    for (int m = 0; m < nl; ++ m)
                    {
                        dq(i, j, k, m) = (dq(i, j, k, m) - rhs0[m]) / diagonal(i, j, k, m);
                    }
                }
            }
        }
        delete[]dq_i;
        delete[]dq_j;
        delete[]dq_k;

        delete[]prim_i;
        delete[]prim_j;
        delete[]prim_k;

        delete[]q_loc_i;
        delete[]q_loc_j;
        delete[]q_loc_k;

        delete[]gykb_i;
        delete[]gykb_j;
        delete[]gykb_k;

        DelPointer2(MG_i);
        DelPointer2(MG_j);
        DelPointer2(MG_k);

        delete[]nxyz;
        delete[]d_q;
        delete[]prim;
        delete[]q_loc;
        delete[]rhs0;

        delete[]de;
        delete[]df;
        delete[]dg;

        delete[]q_pre;
        delete[]dq_pre;
        DelPointer2(MG);
        delete[]temperature;
    }

    void NSSolverStruct::SolveLUSGSBackward(Grid *gridIn, FieldProxy *dq_proxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
    {
        PrintToWindow("Warning: need to check the optimized code, after consider multi-sweep!");
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &dq = dq_proxy->GetField_STR();
        RDouble4D &gxfn = *(grid->GetFaceNormalX());
        RDouble4D &gyfn = *(grid->GetFaceNormalY());
        RDouble4D &gzfn = *(grid->GetFaceNormalZ());
        RDouble4D &gvgn = *(grid->GetFaceNormalVelocity());
        RDouble4D &garea = *(grid->GetFaceArea());
        RDouble3D &gvol = *(grid->GetCellVolume());

        int *iBlank = grid->GetCellTypeContainer();

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
        RDouble refGama = parameters->GetRefGama();
        int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

        RDouble *nxyz = new RDouble[4];
        RDouble *d_q = new RDouble[nEquation];
        RDouble *prim = new RDouble[nEquation];
        RDouble *q_loc = new RDouble[nEquation];
        RDouble *rhs0 = new RDouble[nEquation];

        RDouble *de = new RDouble[nEquation];
        RDouble *df = new RDouble[nEquation];
        RDouble *dg = new RDouble[nEquation];
        RDouble **MG = NewPointer2 <RDouble>(nEquation, nEquation);
        RDouble *temperature = new RDouble[nTemperatureModel];

        RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
        RDouble4D &diagonal = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));

        RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &dt = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));

        RDouble3D &gamma = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));

        Range II, JJ, KK;
        GetRange(ni, nj, nk, -2, 1, II, JJ, KK);
        Range MM(0, nEquation - 1);

        RDouble wmig = 1.0;
        RDouble beta = 1.0;
        RDouble coevis = 0.0;

        Range I, J, K;
        GetRange(ni, nj, nk, 0, -1, I, J, K);

        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        //clz begin
        int viscousType = parameters->GetViscousType();
        //clz end

        //! Obtain the viscosity coefficient, Modified by LiPeng on Mar. 15, 2019.
        RDouble3D viscousLaminar, viscousTurbulence;
        if (viscousType > INVISCID)
        {
            viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
            viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
        }

        //! End of the modification.

        int nDim = GetDim();

        MXDQ_STD_LP(prim, temperature, xfn, yfn, zfn, area, vgn, gama, d_q, df, nNSEquation, nl, rad, -1);
        MXDQ_STD(prim, xfn, yfn, zfn, area, vgn, gama, d_q, df, nNSEquation, nl, rad, -1);
        RDouble *dq_i = new RDouble[nEquation];
        RDouble *dq_j = new RDouble[nEquation];
        RDouble *dq_k = new RDouble[nEquation];

        RDouble *prim_i = new RDouble[nEquation];
        RDouble *prim_j = new RDouble[nEquation];
        RDouble *prim_k = new RDouble[nEquation];

        RDouble *q_loc_i = new RDouble[nEquation];
        RDouble *q_loc_j = new RDouble[nEquation];
        RDouble *q_loc_k = new RDouble[nEquation];

        RDouble *gykb_i = new RDouble[5];
        RDouble *gykb_j = new RDouble[5];
        RDouble *gykb_k = new RDouble[5];

        RDouble **MG_i = NewPointer2 <RDouble>(nEquation, nEquation);
        RDouble **MG_j = NewPointer2 <RDouble>(nEquation, nEquation);
        RDouble **MG_k = NewPointer2 <RDouble>(nEquation, nEquation);


        for (int k = ked; k >= kst; --k)
        {
            int kp1 = k + 1;
            int kp1p1 = kp1 + 1;
            if (nDim != THREE_D)
            {
                kp1 = 1;
                kp1p1 = 1;
            }
            for (int j = jed; j >= jst; --j)
            {
                int jp1 = j + 1;
                int jp1p1 = jp1 + 1;
                for (int i = ied; i >= ist; --i)
                {
                    int ip1 = i + 1;
                    int ip1p1 = ip1 + 1;
                    int cellLabel = (ni - 1) * (nj - 1) * (k - 1) + (ni - 1) * (j - 1) + (i - 1);
                    if (iBlank[cellLabel] == IN_HOLE_COLOR)
                    {
                        for (int m = 0; m < nl; ++ m)
                        {
                            dq(i, j, k, m) = 0.0;
                            //rhs0[m] = 0.0;
                        }
                        continue;
                    }
                    for (int m = 0; m < nl; ++ m)
                    {
                        rhs0[m] = 0.0;
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prim_i[m] = q(ip1, j, k, m);
                        prim_j[m] = q(i, jp1, k, m);
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        q_loc_i[m] = q(ip1, j, k, m);
                        q_loc_j[m] = q(i, jp1, k, m);
                    }

                    for (int m = 0; m < nl; ++ m)
                    {
                        dq_i[m] = dq(ip1, j, k, m);
                        dq_j[m] = dq(i, jp1, k, m);
                    }
                    gykb_i[0] = half * (gxfn(ip1, j, k, 1) + gxfn(ip1p1, j, k, 1));
                    gykb_i[1] = half * (gyfn(ip1, j, k, 1) + gyfn(ip1p1, j, k, 1));
                    gykb_i[2] = 0.0;
                    gykb_i[3] = half * (garea(ip1, j, k, 1) + garea(ip1p1, j, k, 1));
                    gykb_i[4] = half * (gvgn(ip1, j, k, 1) + gvgn(ip1p1, j, k, 1));

                    gykb_j[0] = half * (gxfn(i, jp1, k, 2) + gxfn(i, jp1p1, k, 2));
                    gykb_j[1] = half * (gyfn(i, jp1, k, 2) + gyfn(i, jp1p1, k, 2));
                    gykb_j[2] = 0.0;
                    gykb_j[3] = half * (garea(i, jp1, k, 2) + garea(i, jp1p1, k, 2));
                    gykb_j[4] = half * (gvgn(i, jp1, k, 2) + gvgn(i, jp1p1, k, 2));

                    if (nDim == THREE_D)
                    {
                        gykb_i[2] = half * (gzfn(ip1, j, k, 1) + gzfn(ip1p1, j, k, 1));
                        gykb_j[2] = half * (gzfn(i, jp1, k, 2) + gzfn(i, jp1p1, k, 2));
                    }

                    RDouble &ra = invSpectralRadius(ip1, j, k, 0);
                    RDouble &rb = invSpectralRadius(i, jp1, k, 1);

                    RDouble gama1 = gamma(ip1, j, k);
                    RDouble gama2 = gamma(i, jp1, k);

                    RDouble tm_i = t(ip1, j, k, 0);
                    RDouble tm_j = t(i, jp1, k, 0);
                    if (ifLowSpeedPrecon == 0)
                    {
                        //! The previous method of computing the (M + R) * dQ.
                        MXDQ_STD(prim_i, gykb_i, dq_i, de, gama1, nm, nl, ra, -1);
                        MXDQ_STD(prim_j, gykb_j, dq_j, df, gama2, nm, nl, rb, -1);
                    }
                    else
                    {
                        PreconMatrixFreeDF(q_loc_i, dq_i, gykb_i, MG_i, tm_i, refGama, de, nl, -1);
                        PreconMatrixFreeDF(q_loc_j, dq_j, gykb_j, MG_j, tm_j, refGama, df, nl, -1);
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        dg[m] = 0.0;
                    }

                    if (nDim == THREE_D)
                    {
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            prim_k[m] = q(i, j, kp1, m);
                        }
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            q_loc_k[m] = q(i, j, kp1, m);
                        }
                        for (int m = 0; m < nl; ++ m)
                        {
                            dq_k[m] = dq(i, j, kp1, m);
                        }
                        gykb_k[0] = half * (gxfn(i, j, kp1, 3) + gxfn(i, j, kp1p1, 3));
                        gykb_k[1] = half * (gyfn(i, j, kp1, 3) + gyfn(i, j, kp1p1, 3));
                        gykb_k[2] = half * (gzfn(i, j, kp1, 3) + gzfn(i, j, kp1p1, 3));
                        gykb_k[3] = half * (garea(i, j, kp1, 3) + garea(i, j, kp1p1, 3));
                        gykb_k[4] = half * (gvgn(i, j, kp1, 3) + gvgn(i, j, kp1p1, 3));
                        RDouble &rc = invSpectralRadius(i, j, kp1, 2);
                        RDouble gama3 = gamma(i, j, kp1);
                        RDouble tm_k = t(i, j, kp1, 0);
                        if (ifLowSpeedPrecon == 0)
                        {
                            MXDQ_STD(prim_k, gykb_k, dq_k, dg, gama3, nm, nl, rc, -1);
                        }
                        else
                        {
                            PreconMatrixFreeDF(q_loc_k, dq_k, gykb_k, MG_k, tm_k, refGama, de, nl, -1);
                        }

                    }


                    //! Compute (M - R) * dQ including the real gas conditions, modified by LiPeng on Jan 21, 2019.
                    //MXDQ_STD_LP(prim,xfn,yfn,zfn,area,vgn,gama,d_q,df,nm,nl,rad,-1);
                    //! Solve rhs0.
                    for (int m = 0; m < nl; ++ m)
                    {
                        rhs0[m] += wmig * (de[m] + df[m] + dg[m]);
                    }
                    if (viscousType > INVISCID)    //! Implicit treatment of viscosity term,added by clz : 2012-7-19，modified by LiPeng on Jan 21, 2019.
                    {

                        RDouble rad_vis_i = visSpectralRadius(ip1, j, k, 0);    //! The viscous spectrum radius.
                        RDouble rad_vis_j = visSpectralRadius(i, jp1, k, 1);    //! The viscous spectrum radius.
                        //RDouble rad_vis_k = visSpectralRadius(i  , j  , kp1, 2);    //! The viscous spectrum radius.
                        RDouble rad_vis_k = 0.0;    //! The viscous spectrum radius.

                        if (nDim == THREE_D)
                        {
                            rad_vis_k = visSpectralRadius(i, j, kp1, 2);
                        }

                        for (int m = 0; m < nl; ++ m)
                        {
                            rhs0[m] += -0.5 * (rad_vis_i * dq_i[m] + rad_vis_j * dq_j[m] + rad_vis_k * dq_k[m]);       //! The second method.
                        }
                    }

                    //! Solve dqa.
                    for (int m = 0; m < nl; ++ m)
                    {
                        dq(i, j, k, m) -= rhs0[m] / diagonal(i, j, k, m);
                    }
                }
            }
        }

        delete[] dq_i;
        delete[] dq_j;
        delete[] dq_k;

        delete[] prim_i;
        delete[] prim_j;
        delete[] prim_k;

        delete[] q_loc_i;
        delete[] q_loc_j;
        delete[] q_loc_k;

        delete[] gykb_i;
        delete[] gykb_j;
        delete[] gykb_k;

        DelPointer2(MG_i);
        DelPointer2(MG_j);
        DelPointer2(MG_k);

        //fclose(pfile);

        delete[] nxyz;
        delete[] d_q;
        delete[] prim;
        delete[] q_loc;
        delete[] rhs0;

        delete[] de;
        delete[] df;
        delete[] dg;
        DelPointer2(MG);
        delete[] temperature;
    }
#else
    //! Without any array optimization.
    void NSSolverStruct::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();

        RDouble4D &dq = dqProxy->GetField_STR();
        RDouble4D &deltaFlux = LUplusDQ->GetField_STR();

        RDouble4D &gxfn = *(grid->GetFaceNormalX());
        RDouble4D &gyfn = *(grid->GetFaceNormalY());
        RDouble4D &gzfn = *(grid->GetFaceNormalZ());
        RDouble4D &gvgn = *(grid->GetFaceNormalVelocity());
        RDouble4D &garea = *(grid->GetFaceArea());

        int *iBlank = grid->GetCellTypeContainer();

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
        int nSpeciesNumber = parameters->GetNumberOfSpecies();
        int nEnergyRecycle = parameters->GetnEnergyRecycle();

        int nElectronIndex = parameters->GetIndexOfElectron();

        RDouble refGama = parameters->GetRefGama();

        int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

        RDouble *nxyz = new RDouble[4];
        RDouble *dqNeighbor = new RDouble[nEquation];
        RDouble *prim = new RDouble[nEquation];
        RDouble *qNeighbor = new RDouble[nEquation];
        RDouble *rhs0 = new RDouble[nEquation];
        RDouble *de = new RDouble[nEquation];
        RDouble *df = new RDouble[nEquation];
        RDouble *dg = new RDouble[nEquation];
        RDouble *dqCell = new RDouble[nEquation];
        RDouble *dqOld = new RDouble[nEquation];
        RDouble **MG = NewPointer2 <RDouble>(nEquation, nEquation);
        RDouble *temperature = new RDouble[nTemperatureModel];
        RDouble *qPrecondition = new RDouble[nEquation];
        RDouble *Ux = new RDouble[nEquation];

        RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
        RDouble4D &diagonal = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
        RDouble4D &res = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &gamma = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));

        //RDouble *speciesCvs, RDouble *speciesEtrs, RDouble totalH,  RDouble totalCv, RDouble totalCvtr, RDouble totalCvv, RDouble totalCve
        RDouble4D *speciesCvs = NULL, *speciesEtrs = NULL;
        RDouble3D *totalEnthalpy = NULL, *totalCv = NULL;
        Thermo_Energy *Thermo_Energy_temparay = NULL;
        RDouble3D *totalCvtr = NULL, *totalCvv = NULL, *totalCve = NULL;
        RDouble *Cvs = NULL, *Etrs = NULL;
        if (nEnergyRecycle > 0)
        {
            Thermo_Energy_temparay = gas->GetThermo_Energy_temparay();
            Cvs = Thermo_Energy_temparay->Cvs;
            Etrs = Thermo_Energy_temparay->Etrs;

            speciesCvs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesCvs"));
            speciesEtrs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEtrs"));

            totalEnthalpy = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalEnthalpy"));
            totalCv = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCv"));

            totalCvtr = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCvtr"));
            totalCvv = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCvv"));
            totalCve = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCve"));
        }

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        bool isViscous = parameters->IsViscous();
        RDouble3D viscousLaminar, viscousTurbulence;
        if (isViscous)
        {
            viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
            viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
        }

        int nDim = GetDim();

        for (int k = kCellStart; k <= kCellEnd; ++ k)
        {
            for (int j = jCellStart; j <= jCellEnd; ++ j)
            {
                for (int i = iCellStart; i <= iCellEnd; ++ i)
                {
                    int cellLabel = (ni - 1) * (nj - 1) * (k - 1) + (ni - 1) * (j - 1) + (i - 1);
                    if (iBlank[cellLabel] == IN_HOLE_COLOR)
                    {
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            dq(i, j, k, m) = 0.0;
                        }
                        continue;
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        rhs0[m] = 0.0;

                        //! Temporary variable, only because it is convenient for in/out in function.
                        qPrecondition[m] = q(i, j, k, m);

                        //! Back up the old dq, to compute the convergence.
                        //! it is initialized by zero or res in the first sweep,
                        //! then it is updated by backward sweep, that is dq*.
                        dqOld[m] = dq(i, j, k, m);

                        //! Here, the deltaFlux is computed in Upper backward sweep!
                        //! res: b      deltaFlux: Ux, which is computed in backward sweep.
                        //! the dq is not the real dq, now.
                        //! Although the 'dq' changed can not the right one, it dosen't matter, since 
                        //! the following only using the lower neighbor cells, whose 'dq' has been updated.
                        //! It is convenient for precondition transform.
                        dq(i, j, k, m) = res(i, j, k, m);

                        Ux[m] = deltaFlux(i, j, k, m);

                        //! Then reset it to zero to store the Lower forward sweep!
                        deltaFlux(i, j, k, m) = 0.0;

                        //! Temporary variable, only because it is convenient for in/out in precondition function.
                        dqCell[m] = dq(i, j, k, m);
                    }

                    if (ifLowSpeedPrecon != 0)
                    {
                        RDouble tm = t(i, j, k, IDX::ITT);
                        ComputePreconMatrix(1, MG, qPrecondition, tm, refGama);
                    }

                    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                    {
                        int il1, jl1, kl1, il, jl, kl;
                        grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

                        grid->GetLeftCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, iSurface);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            prim[m] = q(il, jl, kl, m);
                            qNeighbor[m] = q(il, jl, kl, m);
                            dqNeighbor[m] = dq(il, jl, kl, m);
                        }


                        RDouble &xfn1 = gxfn(i, j, k, iSurface);
                        RDouble &yfn1 = gyfn(i, j, k, iSurface);
                        RDouble &zfn1 = gzfn(i, j, k, iSurface);
                        RDouble &area1 = garea(i, j, k, iSurface);
                        RDouble &vgn1 = gvgn(i, j, k, iSurface);

                        RDouble &xfn2 = gxfn(il, jl, kl, iSurface);
                        RDouble &yfn2 = gyfn(il, jl, kl, iSurface);
                        RDouble &zfn2 = gzfn(il, jl, kl, iSurface);
                        RDouble &area2 = garea(il, jl, kl, iSurface);
                        RDouble &vgn2 = gvgn(il, jl, kl, iSurface);

                        RDouble xfn = half * (xfn1 + xfn2);
                        RDouble yfn = half * (yfn1 + yfn2);
                        RDouble zfn = half * (zfn1 + zfn2);
                        RDouble area = half * (area1 + area2);
                        RDouble vgn = half * (vgn1 + vgn2);

                        RDouble &rad = invSpectralRadius(il, jl, kl, iSurface - 1);

                        RDouble gama = gamma(il, jl, kl);

                        if (ifLowSpeedPrecon == 0)
                        {
                            if (nChemical > 0)
                            {
                                for (int m = 0; m < nTemperatureModel; ++ m)
                                {
                                    temperature[m] = t(il, jl, kl, m);
                                }

                                if (nEnergyRecycle == 0)
                                {
                                    ChemicalMXDQ(prim, temperature, xfn, yfn, zfn, area, vgn, dqNeighbor, df, nNSEquation, nLaminar, nTemperatureModel, rad, 1, nElectronIndex);
                                }
                                else
                                {
                                    for (int s = 0; s < nSpeciesNumber; ++ s)
                                    {
                                        Etrs[s] = (*speciesEtrs)(i, j, k, s);
                                        Cvs[s] = (*speciesCvs)(i, j, k, s);
                                    }
                                    ChemicalMXDQR(prim, temperature, xfn, yfn, zfn, area, vgn, dqNeighbor, df, nNSEquation, nLaminar, nTemperatureModel, rad, 1, nElectronIndex,
                                        Cvs, Etrs, (*totalEnthalpy)(i, j, k), (*totalCv)(i, j, k), (*totalCvtr)(i, j, k), (*totalCvv)(i, j, k), (*totalCve)(i, j, k));
                                }
                            }
                            else
                            {
                                MXDQ_STD(prim, xfn, yfn, zfn, area, vgn, gama, dqNeighbor, df, nNSEquation, nLaminar, rad, 1);
                            }
                        }
                        else
                        {
                            PreconMatrixFreeDF(qNeighbor, dqNeighbor, xfn, yfn, zfn, area, vgn, MG, refGama, df, nLaminar, 1);
                        }

                        //! Add Flux together.
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            rhs0[m] -= df[m];
                        }

                        if (isViscous)
                        {
                            RDouble rad_vis = half * visSpectralRadius(il, jl, kl, iSurface - 1);

                            for (int m = 0; m < nEquation; ++ m)
                            {
                                rhs0[m] -= rad_vis * dqNeighbor[m];
                            }
                        }
                    }

                    if (ifLowSpeedPrecon == 1)
                    {
                        //! Here, the 'dqCell' is 'res' actually.
                        //! The reason why operate on 'dqCell' other than 'res', is that the 'dqCell' is a temporary variable.
                        //! and the real 'res' is forbidden to be changed.
                        PreconMatrixDF(dqCell, MG, nLaminar);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            dq(i, j, k, m) = dqCell[m];
                        }
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        //! dq = { (b - Ux) - Lx } / D.
                        //! Note: the 'dq' has actually initialized by the residual at the beginning of this function.
                        //! the 'rhs0' is the total Delta Flux of the neighbors.            
                        //! rhs0: Lx    diagonal: D.
                        dq(i, j, k, m) = (dq(i, j, k, m) - Ux[m] - rhs0[m]) / diagonal(i, j, k, m);
                        //! Store the lower forward sweep delta-flux, which will be used in the backward sweep.
                        deltaFlux(i, j, k, m) += rhs0[m];

                        sweepNormal += SQR(dq(i, j, k, m) - dqOld[m]);
                    }

                }
            }
        }

        delete[] nxyz;    nxyz = NULL;
        delete[] dqNeighbor;    dqNeighbor = NULL;
        delete[] prim;    prim = NULL;
        delete[] qNeighbor;    qNeighbor = NULL;
        delete[] rhs0;    rhs0 = NULL;

        delete[] de;    de = NULL;
        delete[] df;    df = NULL;
        delete[] dg;    dg = NULL;

        delete[] qPrecondition;    qPrecondition = NULL;
        delete[] dqCell;    dqCell = NULL;
        DelPointer2(MG);
        delete[] temperature;    temperature = NULL;

        delete[] dqOld;    dqOld = NULL;
        delete[] Ux;    Ux = NULL;
    }

    void NSSolverStruct::SolveLUSGSBackward(Grid *gridIn, FieldProxy *dq_proxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        RDouble4D &dq = dq_proxy->GetField_STR();
        RDouble4D &deltaFlux = LUplusDQ->GetField_STR();

        RDouble4D &gxfn = *(grid->GetFaceNormalX());
        RDouble4D &gyfn = *(grid->GetFaceNormalY());
        RDouble4D &gzfn = *(grid->GetFaceNormalZ());
        RDouble4D &gvgn = *(grid->GetFaceNormalVelocity());
        RDouble4D &garea = *(grid->GetFaceArea());

        int *iBlank = grid->GetCellTypeContainer();

        Param_NSSolverStruct *parameters = GetControlParameters();
        int nNSEquation = parameters->GetNSEquationNumber();
        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
        int nSpeciesNumber = parameters->GetNumberOfSpecies();
        int nEnergyRecycle = parameters->GetnEnergyRecycle();

        int nElectronIndex = parameters->GetIndexOfElectron();

        RDouble refGama = parameters->GetRefGama();
        int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

        RDouble *nxyz = new RDouble[4];
        RDouble *dqNeighbor = new RDouble[nEquation];
        RDouble *prim = new RDouble[nEquation];
        RDouble *qNeighbor = new RDouble[nEquation];
        RDouble *dqOld = new RDouble[nEquation];
        RDouble *rhs0 = new RDouble[nEquation];
        RDouble *de = new RDouble[nEquation];
        RDouble *df = new RDouble[nEquation];
        RDouble *dg = new RDouble[nEquation];
        RDouble **MG = NewPointer2 <RDouble>(nEquation, nEquation);
        RDouble *temperature = new RDouble[nTemperatureModel];
        RDouble *qPrecondition = new RDouble[nEquation];
        RDouble *dqCell = new RDouble[nEquation];
        RDouble *Lx = new RDouble[nEquation];

        RDouble4D &invSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
        RDouble4D &diagonal = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
        RDouble4D &res = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &visSpectralRadius = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &gamma = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));

        //RDouble *speciesCvs, RDouble *speciesEtrs, RDouble totalH,  RDouble totalCv, RDouble totalCvtr, RDouble totalCvv, RDouble totalCve
        RDouble4D *speciesCvs = NULL, *speciesEtrs = NULL;
        RDouble3D *totalEnthalpy = NULL, *totalCv = NULL;
        Thermo_Energy *Thermo_Energy_temparay = NULL;
        RDouble3D *totalCvtr = NULL, *totalCvv = NULL, *totalCve = NULL;
        RDouble *Cvs = NULL, *Etrs = NULL;
        if (nEnergyRecycle > 0)
        {
            Thermo_Energy_temparay = gas->GetThermo_Energy_temparay();
            Cvs = Thermo_Energy_temparay->Cvs;
            Etrs = Thermo_Energy_temparay->Etrs;

            speciesCvs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesCvs"));
            speciesEtrs = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("speciesEtrs"));

            totalEnthalpy = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalEnthalpy"));
            totalCv = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCv"));

            totalCvtr = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCvtr"));
            totalCvv = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCvv"));
            totalCve = reinterpret_cast <RDouble3D *> (grid->GetDataPtr("totalCve"));
        }

        Range II, JJ, KK;
        GetRange(ni, nj, nk, -2, 1, II, JJ, KK);
        Range MM(0, nEquation - 1);

        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        grid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, NON_GHOST);

        bool isViscous = parameters->IsViscous();
        RDouble3D viscousLaminar, viscousTurbulence;
        if (isViscous)
        {
            viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
            viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
        }

        int nDim = GetDim();

        for (int k = kCellEnd; k >= kCellStart; --k)
        {
            for (int j = jCellEnd; j >= jCellStart; --j)
            {
                for (int i = iCellEnd; i >= iCellStart; --i)
                {
                    int cellLabel = (ni - 1) * (nj - 1) * (k - 1) + (ni - 1) * (j - 1) + (i - 1);
                    if (iBlank[cellLabel] == IN_HOLE_COLOR)
                    {
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            dq(i, j, k, m) = 0.0;
                        }
                        continue;
                    }
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        rhs0[m] = 0.0;
                        dqNeighbor[m] = dq(i, j, k, m);

                        //! Temporary variable, only because it is convenient for in/out in function.
                        qPrecondition[m] = q(i, j, k, m);

                        //! Back up the old dq, to compute the convergence.
                        //! the 'dq' is dq*, which has been updated in forward.
                        dqOld[m] = dq(i, j, k, m);

                        //! the dq is not the real dq, now.
                        //! it is convenient for precondition transform.
                        dq(i, j, k, m) = res(i, j, k, m);

                        Lx[m] = deltaFlux(i, j, k, m);

                        //! Then reset it to zero to store the upper forward sweep!
                        deltaFlux(i, j, k, m) = 0.0;

                        //! Temporary variable, only because it is convenient for in/out in precondition function.
                        dqCell[m] = dq(i, j, k, m);
                    }

                    if (ifLowSpeedPrecon != 0)
                    {
                        RDouble tm = t(i, j, k, IDX::ITT);
                        ComputePreconMatrix(1, MG, qPrecondition, tm, refGama);
                    }

                    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                    {
                        int il1, jl1, kl1, il, jl, kl;
                        grid->GetNsurfIndex(il1, jl1, kl1, iSurface);

                        grid->GetRightCellOfFace(i, j, k, il1, jl1, kl1, il, jl, kl, iSurface);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            prim[m] = q(il, jl, kl, m);
                            qNeighbor[m] = q(il, jl, kl, m);
                            dqNeighbor[m] = dq(il, jl, kl, m);
                        }

                        //clz begin
                        RDouble &xfn1 = gxfn(il, jl, kl, iSurface);
                        RDouble &yfn1 = gyfn(il, jl, kl, iSurface);
                        RDouble &zfn1 = gzfn(il, jl, kl, iSurface);
                        RDouble &area1 = garea(il, jl, kl, iSurface);
                        RDouble &vgn1 = gvgn(il, jl, kl, iSurface);

                        RDouble &xfn2 = gxfn(il + il1, jl + jl1, kl + kl1, iSurface);
                        RDouble &yfn2 = gyfn(il + il1, jl + jl1, kl + kl1, iSurface);
                        RDouble &zfn2 = gzfn(il + il1, jl + jl1, kl + kl1, iSurface);
                        RDouble &area2 = garea(il + il1, jl + jl1, kl + kl1, iSurface);
                        RDouble &vgn2 = gvgn(il + il1, jl + jl1, kl + kl1, iSurface);

                        RDouble xfn = half * (xfn1 + xfn2);
                        RDouble yfn = half * (yfn1 + yfn2);
                        RDouble zfn = half * (zfn1 + zfn2);
                        RDouble area = half * (area1 + area2);
                        RDouble vgn = half * (vgn1 + vgn2);

                        RDouble &rad = invSpectralRadius(il, jl, kl, iSurface - 1);

                        RDouble gama = gamma(il, jl, kl);

                        if (ifLowSpeedPrecon == 0)
                        {
                            //! Compute (M - R) * dQ including the real gas conditions.
                            if (nChemical > 0)
                            {
                                for (int m = 0; m < nTemperatureModel; ++ m)
                                {
                                    temperature[m] = t(il, jl, kl, m);
                                }
                                if (nEnergyRecycle == 0)
                                {
                                    ChemicalMXDQ(prim, temperature, xfn, yfn, zfn, area, vgn, dqNeighbor, df, nNSEquation, nLaminar, nTemperatureModel, rad, -1, nElectronIndex);
                                }
                                else
                                {
                                    for (int s = 0; s < nSpeciesNumber; ++ s)
                                    {
                                        Etrs[s] = (*speciesEtrs)(i, j, k, s);
                                        Cvs[s] = (*speciesCvs)(i, j, k, s);
                                    }
                                    ChemicalMXDQR(prim, temperature, xfn, yfn, zfn, area, vgn, dqNeighbor, df, nNSEquation, nLaminar, nTemperatureModel, rad, -1, nElectronIndex,
                                        Cvs, Etrs, (*totalEnthalpy)(i, j, k), (*totalCv)(i, j, k), (*totalCvtr)(i, j, k), (*totalCvv)(i, j, k), (*totalCve)(i, j, k));
                                }
                            }
                            else
                            {
                                MXDQ_STD(prim, xfn, yfn, zfn, area, vgn, gama, dqNeighbor, df, nNSEquation, nLaminar, rad, -1);
                            }
                        }
                        else
                        {
                            PreconMatrixFreeDF(qNeighbor, dqNeighbor, xfn, yfn, zfn, area, vgn, MG, refGama, df, nLaminar, -1);
                        }

                        //! Add Flux together.
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            rhs0[m] += df[m];
                        }

                        if (isViscous)
                        {
                            RDouble rad_vis = half * visSpectralRadius(il, jl, kl, iSurface - 1);

                            for (int m = 0; m < nEquation; ++ m)
                            {
                                rhs0[m] += -rad_vis * dqNeighbor[m];
                            }
                        }
                    }

                    if (ifLowSpeedPrecon == 1)
                    {
                        //! Here, the 'dqCell' is 'res' actually.
                        //! The reason why operate on 'dqCell' other than 'res', is that the 'dqCell' is a temporary variable.
                        //! and the real 'res' is forbidden to be changed.
                        PreconMatrixDF(dqCell, MG, nLaminar);

                        for (int m = 0; m < nEquation; ++ m)
                        {
                            dq(i, j, k, m) = dqCell[m];
                        }
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        //! Note: the 'dq' has been updated by the forward sweep.
                        //! the 'rhs0' is the total Delta Flux of the neighbors.
                        //! x = {(b - LX) - Ux} / D.
                        //! rhs0: Ux.    diagonal: D.
                        dq(i, j, k, m) = (dq(i, j, k, m) - Lx[m] - rhs0[m]) / diagonal(i, j, k, m);

                        //! Store the upper backward sweep delta-flux, which will be used in the forward sweep.
                        deltaFlux(i, j, k, m) += rhs0[m];

                        sweepNormal += SQR(dq(i, j, k, m) - dqOld[m]);
                    }
                }
            }
        }

        delete[] nxyz;    nxyz = NULL;
        delete[] dqNeighbor;    dqNeighbor = NULL;
        delete[] dqOld;    dqOld = NULL;
        delete[] prim;    prim = NULL;
        delete[] qNeighbor;    qNeighbor = NULL;
        delete[] rhs0;    rhs0 = NULL;

        delete[] de;    de = NULL;
        delete[] df;    df = NULL;
        delete[] dg;    dg = NULL;
        DelPointer2(MG);
        delete[] temperature;    temperature = NULL;
        delete[] Lx;    Lx = NULL;
        delete[] qPrecondition;    qPrecondition = NULL;
        delete[] dqCell;    dqCell = NULL;
    }
#endif

#ifdef USE_INHOUSE_VERSION
    void NSSolverStruct::SolveBLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
    {
        StructGrid *grid = StructGridCast(gridIn);
        SolveBLUSGSForwardInHouse(grid, dqProxy, LUplusDQ, sweepNormal);
    }

    void NSSolverStruct::SolveBLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
    {
        StructGrid *grid = StructGridCast(gridIn);
        SolveBLUSGSBackwardInHouse(grid, dqProxy, LUplusDQ, sweepNormal);

    }
#endif

    void NSSolverStruct::SetGhostDQLUSGS(Grid *gridIn, RDouble4D &dq)
    {
        StructGrid *grid = StructGridCast(gridIn);

        RDouble4D &xfn = *(grid->GetFaceNormalX());
        RDouble4D &yfn = *(grid->GetFaceNormalY());
        RDouble4D &zfn = *(grid->GetFaceNormalZ());

        Param_NSSolverStruct *parameters = GetControlParameters();

        int nLaminar = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

        int mTV = parameters->GetmTV();
        int mTE = parameters->GetmTE();

        int viscousType = parameters->GetViscousType();

        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &t = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble3D &gama = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));

        using namespace IDX;
        using namespace GAS_SPACE;
        using namespace PHENGLEI;
        RDouble *prim = new RDouble[nEquation];
        RDouble *Q = new RDouble[nEquation];

        CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();

            int *s_lr3d = structBC->GetFaceDirectionIndex();
            int iSurface = structBC->GetFaceDirection() + 1;

            int ist, ied, jst, jed, kst, ked;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        int is = i;
                        int js = j;
                        int ks = k;

                        int it = i + s_lr3d[0];
                        int jt = j + s_lr3d[1];
                        int kt = k + s_lr3d[2];

                        //! Conservative Variables on level n+1.
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            prim[m] = q(is, js, ks, m);
                        }

                        RDouble tv = t(is, js, ks, mTV);
                        RDouble te = t(is, js, ks, mTE);

                        gas->Primitive2Conservative(prim, gama(is, js, ks), tv, te, Q);

                        for (int m = 0; m < nLaminar; ++ m)
                        {
                            Q[m] += dq(is, js, ks, m);
                        }

                        //! Non-conservative variables in re in level n.
                        RDouble rm = q(it, jt, kt, IR);
                        RDouble um = q(it, jt, kt, IU);
                        RDouble vm = q(it, jt, kt, IV);
                        RDouble wm = q(it, jt, kt, IW);
                        RDouble pm = q(it, jt, kt, IP);

                        int in, jn, kn;
                        structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                        RDouble nxs = xfn(in, jn, kn, iSurface);
                        RDouble nys = yfn(in, jn, kn, iSurface);
                        RDouble nzs = zfn(in, jn, kn, iSurface);

                        if (BCType == PHENGLEI::SYMMETRY || (IsWall(BCType) && viscousType == INVISCID))
                        {
                            //! Inviscid wall or symmetry boundary.
                            dq(it, jt, kt, IR) = dq(is, js, ks, IR);
                            RDouble rvn = two * (nxs * dq(is, js, ks, IU) + nys * dq(is, js, ks, IV) + nzs * dq(is, js, ks, IW));
                            dq(it, jt, kt, IU) = dq(is, js, ks, IU) - nxs * rvn;
                            dq(it, jt, kt, IV) = dq(is, js, ks, IV) - nys * rvn;
                            dq(it, jt, kt, IW) = dq(is, js, ks, IW) - nzs * rvn;
                            rvn = two * (nxs * Q[IU] + nys * Q[IV] + nzs * Q[IW]);
                            RDouble ru = Q[IU] - nxs * rvn;
                            RDouble rv = Q[IV] - nys * rvn;
                            RDouble rw = Q[IW] - nzs * rvn;
                            RDouble p_new = Q[IP] - half * (Q[IU] * Q[IU] + Q[IV] * Q[IV] + Q[IW] * Q[IW]) / Q[IR];

                            RDouble gam1 = gama(is, js, ks) - one;
                            dq(it, jt, kt, IP) = p_new - pm / gam1 + half * ((ru * ru + rv * rv + rw * rw) / Q[IR] - rm * (um * um + vm * vm + wm * wm));
                        }
                        else if (IsWall(BCType) && viscousType > INVISCID)
                        {
                            RDouble ub = 0.0;
                            RDouble vb = 0.0;
                            RDouble wb = 0.0;
                            dq(it, jt, kt, IR) = dq(is, js, ks, IR);
                            dq(it, jt, kt, IU) = -dq(is, js, ks, IU) + two * ub * dq(is, js, ks, IR);
                            dq(it, jt, kt, IV) = -dq(is, js, ks, IV) + two * vb * dq(is, js, ks, IR);
                            dq(it, jt, kt, IW) = -dq(is, js, ks, IW) + two * wb * dq(is, js, ks, IR);
                            dq(it, jt, kt, IP) = dq(is, js, ks, IP) + ub * (dq(it, jt, kt, IU) - dq(is, js, ks, IU))
                                + vb * (dq(it, jt, kt, IV) - dq(is, js, ks, IV))
                                + wb * (dq(it, jt, kt, IW) - dq(is, js, ks, IW));
                        }
                        else if (BCType == PHENGLEI::EXTRAPOLATION)
                        {
                            for (int m = 0; m < nLaminar; ++ m)
                            {
                                dq(it, jt, kt, m) = dq(is, js, ks, m);
                            }
                        }
                        else
                        {
                            //! FIX_ALL set DQ = 0.
                            for (int m = 0; m < nLaminar; ++ m)
                            {
                                dq(it, jt, kt, m) = zero;
                            }
                        }
                    }
                }
            }
        }
        delete[]prim;
        delete[]Q;
    }

    void NSSolverStruct::UploadInterfaceData(ActionKey *actkey)
    {
        Param_NSSolverStruct *parameters = GetControlParameters();
        int level = actkey->level;
        StructGrid *grid = reinterpret_cast <StructGrid *> (GetGrid(level));
        StructGrid *fineGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *interfaceInfo = fineGrid->GetInterfaceInfo();
        if (!interfaceInfo)
        {
            return;
        }

        int nEquation = GetNumberOfEquations();
        RDouble4D *primitiveVariables = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        PHSPACE::UploadInterfaceValue(grid, primitiveVariables, "q", nEquation);

        int nTemperatureModel = parameters->GetTemperatureModel();
        RDouble4D *temperatures = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        PHSPACE::UploadInterfaceValue(grid, temperatures, "t", nTemperatureModel);
    }

    void NSSolverStruct::DownloadInterfaceData(ActionKey *actkey)
    {
        Param_NSSolverStruct *parameters = GetControlParameters();
        int level = actkey->level;
        StructGrid *grid = reinterpret_cast<StructGrid *>(GetGrid(level));
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        if (!interfaceInfo)
        {
            return;
        }

        int nEquation = GetNumberOfEquations();

        RDouble4D *primitiveVariables = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        PHSPACE::DownloadInterfaceValue(grid, primitiveVariables, "q", nEquation);

        int nTemperatureModel = parameters->GetTemperatureModel();
        RDouble4D *temperatures = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        PHSPACE::DownloadInterfaceValue(grid, temperatures, "t", nTemperatureModel);
    }

    void NSSolverStruct::UploadOversetData(ActionKey *actkey)
    {
        int level = actkey->level;
        StructGrid *grid = StructGridCast(GetGrid(level));
        Param_NSSolverStruct *parameters = GetControlParameters();
        int nEquation = GetNumberOfEquations();

        RDouble4D *primitiveVariables = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        PHSPACE::UploadOversetValue(grid, primitiveVariables, "q", nEquation);

        int nTemperatureModel = parameters->GetTemperatureModel();
        RDouble4D *temperatures = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        PHSPACE::UploadOversetValue(grid, temperatures, "t", nTemperatureModel);
    }

    void NSSolverStruct::DownloadOversetData(ActionKey *actkey)
    {
        int level = actkey->level;
        StructGrid *grid = reinterpret_cast<StructGrid *>(GetGrid(level));
        Param_NSSolverStruct *parameters = GetControlParameters();

        int nEquation = GetNumberOfEquations();

        RDouble4D *primitiveVariables = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        PHSPACE::DownloadOversetValue(grid, primitiveVariables, "q", nEquation);

        int nTemperatureModel = parameters->GetTemperatureModel();
        RDouble4D *temperatures = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        PHSPACE::DownloadOversetValue(grid, temperatures, "t", nTemperatureModel);
    }

    void NSSolverStruct::CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int level = grid->GetLevel();
        InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation)
        {
            return;
        }

        StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

        int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
        int *interfaceIndexContainerForSend = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

        dataContainer->moveto_begin();

        if (0 == nEquation)
        {
            RDouble3D &fieldSend = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr(fieldName));

            for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
            {
                for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
                {
                    int iFace = interfaceIndexContainerForSend[iLocalFace];
                    int is, js, ks;
                    finestGrid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
                    finestGrid->RemapMultigridIJK(level, is, js, ks);
                    PHWrite(dataContainer, fieldSend(is, js, ks));
                }
            }
        }
        else
        {
            RDouble4D &fieldSend = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr(fieldName));

            for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
            {
                for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
                {
                    int iFace = interfaceIndexContainerForSend[iLocalFace];
                    int is, js, ks;
                    finestGrid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
                    finestGrid->RemapMultigridIJK(level, is, js, ks);
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        PHWrite(dataContainer, fieldSend(is, js, ks, m));
                    }
                }
            }
        }
    }

    void NSSolverStruct::DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int level = grid->GetLevel();
        InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation)
        {
            return;
        }

        StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

        int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
        int *interfaceIndexContainerForReceive = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

        dataContainer->moveto_begin();

        if (0 == nEquation)
        {
            RDouble3D &fieldRecv = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr(fieldName));

            for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
            {
                for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
                {
                    int iFace = interfaceIndexContainerForReceive[iLocalFace];
                    int it, jt, kt;
                    finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
                    finestGrid->RemapMultigridIJK(level, it, jt, kt);
                    PHRead(dataContainer, fieldRecv(it, jt, kt));
                }
            }
        }
        else
        {
            RDouble4D &fieldRecv = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr(fieldName));

            for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
            {
                for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
                {
                    int iFace = interfaceIndexContainerForReceive[iLocalFace];
                    int it, jt, kt;
                    finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
                    finestGrid->RemapMultigridIJK(level, it, jt, kt);
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        PHRead(dataContainer, fieldRecv(it, jt, kt, m));
                    }
                }
            }
        }
    }

    void NSSolverStruct::RotateVectorFromInterface(Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int level = grid->GetLevel();
        InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation)
        {
            return;
        }

        StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();
        CompositeStructBC *composite_bcregion = finestGrid->GetCompositeBCRegion();
        RDouble4D &rotNSgradValueX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotNSgradValueX"));
        RDouble4D &rotNSgradValueY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotNSgradValueY"));
        RDouble4D &rotNSgradValueZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotNSgradValueZ"));

        rotNSgradValueX = 0.0;
        rotNSgradValueY = 0.0;
        rotNSgradValueZ = 0.0;

        int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
        int *interfaceIndexContainerForReceive = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);
        double rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
        rotationAngle = rotationAngle * PI / 180;
        if (nEquation > 0)
        {
            RDouble4D &fieldRecvNSX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
            RDouble4D &fieldRecvNSY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
            RDouble4D &fieldRecvNSZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

            for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
            {
                for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
                {
                    int iFace = interfaceIndexContainerForReceive[iLocalFace];
                    int it, jt, kt;
                    finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
                    finestGrid->RemapMultigridIJK(level, it, jt, kt);
                    int *ibcregions = composite_bcregion->GetIFaceInfo();
                    int iBCRegion = ibcregions[iFace];
                    StructBC *bcregion = composite_bcregion->GetBCRegion(iBCRegion);
                    string bcName = bcregion->GetBCName();
                    if (bcName == "Periodic_up")
                    {
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            rotNSgradValueY(it, jt, kt, m) = fieldRecvNSY(it, jt, kt, m) * cos(2 * PI - rotationAngle) - fieldRecvNSZ(it, jt, kt, m) * sin(2 * PI - rotationAngle);
                            rotNSgradValueZ(it, jt, kt, m) = fieldRecvNSY(it, jt, kt, m) * sin(2 * PI - rotationAngle) + fieldRecvNSZ(it, jt, kt, m) * cos(2 * PI - rotationAngle);
                        }
                    }
                    else if (bcName == "Periodic_down")
                    {
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            rotNSgradValueY(it, jt, kt, m) = fieldRecvNSY(it, jt, kt, m) * cos(rotationAngle) - fieldRecvNSZ(it, jt, kt, m) * sin(rotationAngle);
                            rotNSgradValueZ(it, jt, kt, m) = fieldRecvNSY(it, jt, kt, m) * sin(rotationAngle) + fieldRecvNSZ(it, jt, kt, m) * cos(rotationAngle);
                        }
                    }

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        fieldRecvNSY(it, jt, kt, m) = rotNSgradValueY(it, jt, kt, m);
                        fieldRecvNSZ(it, jt, kt, m) = rotNSgradValueZ(it, jt, kt, m);
                    }
                }
            }
        }
    }

    void NSSolverStruct::CompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int level = grid->GetLevel();
        InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation)
        {
            return;
        }

        StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

        int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
        int *interfaceIndexContainerForSend = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

        dataContainer->moveto_begin();

        RDouble4D &dqdx_cc = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdx_cc_ruvwpt"));
        RDouble4D &dqdy_cc = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdy_cc_ruvwpt"));
        RDouble4D &dqdz_cc = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("dqdz_cc_ruvwpt"));

        RDouble limitcoefofinterface = GlobalDataBase::GetDoubleParaFromDB("limitcoefofinterface");
        int limitVector = GlobalDataBase::GetIntParaFromDB("limitVector");

        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForSend[iLocalFace];
            int is, js, ks;
            finestGrid->GetSourceIndexIJK(iFace, 1, is, js, ks);
            finestGrid->RemapMultigridIJK(level, is, js, ks);

            for (int m = 0; m <= 5; ++ m)
            {
                PHWrite(dataContainer, dqdx_cc(is, js, ks, m));
                PHWrite(dataContainer, dqdy_cc(is, js, ks, m));
                PHWrite(dataContainer, dqdz_cc(is, js, ks, m));
            }

            if (limitVector == 0)
            {
                //! scalar
                PHWrite(dataContainer, limitcoefofinterface);
            }
            else
            {
                PHWrite(dataContainer, limitcoefofinterface);
                PHWrite(dataContainer, limitcoefofinterface);
                PHWrite(dataContainer, limitcoefofinterface);
                PHWrite(dataContainer, limitcoefofinterface);
                PHWrite(dataContainer, limitcoefofinterface);
            }
        }
    }

    void NSSolverStruct::DecompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int level = grid->GetLevel();
        InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation)
        {
            return;
        }

        StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

        int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
        int *interfaceIndexContainerForReceive = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

        dataContainer->moveto_begin();

        RDouble4D &gradUVWTCellCenterX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
        RDouble4D &gradUVWTCellCenterY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
        RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

        RDouble dqdx_cc_ruvwpt[6];
        RDouble dqdy_cc_ruvwpt[6];
        RDouble dqdz_cc_ruvwpt[6];

        int limitVector = GlobalDataBase::GetIntParaFromDB("limitVector");
        RDouble limitofunstr_nouse = 0.0;

        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForReceive[iLocalFace];
            int it, jt, kt;
            finestGrid->GetTargetIndexIJK(iFace, 1, it, jt, kt);
            finestGrid->RemapMultigridIJK(level, it, jt, kt);
            for (int m = 0; m <= 5; ++ m)
            {
                PHRead(dataContainer, dqdx_cc_ruvwpt[m]);
                PHRead(dataContainer, dqdy_cc_ruvwpt[m]);
                PHRead(dataContainer, dqdz_cc_ruvwpt[m]);
            }

            if (limitVector == 0)
            {
                //! scalar
                PHRead(dataContainer, limitofunstr_nouse);
            }
            else
            {
                PHRead(dataContainer, limitofunstr_nouse);
                PHRead(dataContainer, limitofunstr_nouse);
                PHRead(dataContainer, limitofunstr_nouse);
                PHRead(dataContainer, limitofunstr_nouse);
                PHRead(dataContainer, limitofunstr_nouse);
            }

            gradUVWTCellCenterX(it, jt, kt, 0) = dqdx_cc_ruvwpt[1]; //grad of U
            gradUVWTCellCenterY(it, jt, kt, 0) = dqdy_cc_ruvwpt[1]; //grad of U
            gradUVWTCellCenterZ(it, jt, kt, 0) = dqdz_cc_ruvwpt[1]; //grad of U

            gradUVWTCellCenterX(it, jt, kt, 1) = dqdx_cc_ruvwpt[2]; //grad of V
            gradUVWTCellCenterY(it, jt, kt, 1) = dqdy_cc_ruvwpt[2]; //grad of V
            gradUVWTCellCenterZ(it, jt, kt, 1) = dqdz_cc_ruvwpt[2]; //grad of V

            gradUVWTCellCenterX(it, jt, kt, 2) = dqdx_cc_ruvwpt[3]; //grad of W
            gradUVWTCellCenterY(it, jt, kt, 2) = dqdy_cc_ruvwpt[3]; //grad of W
            gradUVWTCellCenterZ(it, jt, kt, 2) = dqdz_cc_ruvwpt[3]; //grad of W

            gradUVWTCellCenterX(it, jt, kt, 3) = dqdx_cc_ruvwpt[5]; //grad of T
            gradUVWTCellCenterY(it, jt, kt, 3) = dqdy_cc_ruvwpt[5]; //grad of T
            gradUVWTCellCenterZ(it, jt, kt, 3) = dqdz_cc_ruvwpt[5]; //grad of T
        }
    }

    void NSSolverStruct::CompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int level = grid->GetLevel();
        InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation)
        {
            return;
        }

        StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

        int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
        int *interfaceIndexContainerForSend = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

        dataContainer->moveto_begin();

        RDouble4D &QLI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1"));
        RDouble4D &QLJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2"));
        RDouble4D &QLK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3"));
        RDouble4D &QRI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1"));
        RDouble4D &QRJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2"));
        RDouble4D &QRK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3"));

        RDouble qface[5];
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForSend[iLocalFace];

            int is, js, ks;
            finestGrid->GetSourceIndexIJK(iFace, 1, is, js, ks);
            finestGrid->RemapMultigridIJK(level, is, js, ks);

            int id, jd, kd;
            int i_lr, j_lr, k_lr;
            finestGrid->GetSourceIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);

            if (i_lr != 0)
            {
                //! I_face
                if (i_lr == -1)
                {
                    //! left boundary
                    qface[0] = QRI(is, js, ks, 0);
                    qface[1] = QRI(is, js, ks, 1);
                    qface[2] = QRI(is, js, ks, 2);
                    qface[3] = QRI(is, js, ks, 3);
                    qface[4] = QRI(is, js, ks, 4);
                }
                else
                {
                    //! right boundary
                    qface[0] = QLI(is, js, ks, 0);
                    qface[1] = QLI(is, js, ks, 1);
                    qface[2] = QLI(is, js, ks, 2);
                    qface[3] = QLI(is, js, ks, 3);
                    qface[4] = QLI(is, js, ks, 4);
                }
            }
            else if (j_lr != 0)
            {
                //! J_face
                if (j_lr == -1)
                {
                    //! left boundary
                    qface[0] = QRJ(is, js, ks, 0);
                    qface[1] = QRJ(is, js, ks, 1);
                    qface[2] = QRJ(is, js, ks, 2);
                    qface[3] = QRJ(is, js, ks, 3);
                    qface[4] = QRJ(is, js, ks, 4);
                }
                else
                {
                    //! right boundary
                    qface[0] = QLJ(is, js, ks, 0);
                    qface[1] = QLJ(is, js, ks, 1);
                    qface[2] = QLJ(is, js, ks, 2);
                    qface[3] = QLJ(is, js, ks, 3);
                    qface[4] = QLJ(is, js, ks, 4);
                }
            }
            else
            {
                //! K_face
                if (k_lr == -1)
                {
                    //! left boundary
                    qface[0] = QRK(is, js, ks, 0);
                    qface[1] = QRK(is, js, ks, 1);
                    qface[2] = QRK(is, js, ks, 2);
                    qface[3] = QRK(is, js, ks, 3);
                    qface[4] = QRK(is, js, ks, 4);
                }
                else
                {
                    //! right boundary
                    qface[0] = QLK(is, js, ks, 0);
                    qface[1] = QLK(is, js, ks, 1);
                    qface[2] = QLK(is, js, ks, 2);
                    qface[3] = QLK(is, js, ks, 3);
                    qface[4] = QLK(is, js, ks, 4);
                }
            }
            for (int m = 0; m < 5; ++ m)
            {
                PHWrite(dataContainer, qface[m]);
            }
        }
    }

    void NSSolverStruct::DecompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int level = grid->GetLevel();
        InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation)
        {
            return;
        }

        StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

        int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
        int *interfaceIndexContainerForReceive = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

        dataContainer->moveto_begin();

        RDouble4D &QLI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1"));
        RDouble4D &QLJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2"));
        RDouble4D &QLK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3"));
        RDouble4D &QRI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1"));
        RDouble4D &QRJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2"));
        RDouble4D &QRK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3"));

        RDouble qface[5];
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForReceive[iLocalFace];
            int it, jt, kt;
            finestGrid->GetTargetIndexIJK(iFace, 1, it, jt, kt);
            finestGrid->RemapMultigridIJK(level, it, jt, kt);

            for (int m = 0; m < 5; ++ m)
            {
                PHRead(dataContainer, qface[m]);
            }

            int id, jd, kd;
            int i_lr, j_lr, k_lr;
            finestGrid->GetTargetIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);
            if (i_lr != 0)
            {
                //! I_face
                if (i_lr == -1)
                {
                    //! left boundary
                    QLI(it, jt, kt, 0) = qface[0];
                    QLI(it, jt, kt, 1) = qface[1];
                    QLI(it, jt, kt, 2) = qface[2];
                    QLI(it, jt, kt, 3) = qface[3];
                    QLI(it, jt, kt, 4) = qface[4];
                }
                else
                {
                    //! right boundary
                    QRI(it, jt, kt, 0) = qface[0];
                    QRI(it, jt, kt, 1) = qface[1];
                    QRI(it, jt, kt, 2) = qface[2];
                    QRI(it, jt, kt, 3) = qface[3];
                    QRI(it, jt, kt, 4) = qface[4];
                }
            }
            else if (j_lr != 0)
            {
                //! J_face
                if (j_lr == -1)
                {
                    //! left boundary
                    QLJ(it, jt, kt, 0) = qface[0];
                    QLJ(it, jt, kt, 1) = qface[1];
                    QLJ(it, jt, kt, 2) = qface[2];
                    QLJ(it, jt, kt, 3) = qface[3];
                    QLJ(it, jt, kt, 4) = qface[4];
                }
                else
                {
                    //! right boundary
                    QRJ(it, jt, kt, 0) = qface[0];
                    QRJ(it, jt, kt, 1) = qface[1];
                    QRJ(it, jt, kt, 2) = qface[2];
                    QRJ(it, jt, kt, 3) = qface[3];
                    QRJ(it, jt, kt, 4) = qface[4];
                }
            }
            else
            {
                //! K_face
                if (k_lr == -1)
                {
                    //! left boundary
                    QLK(it, jt, kt, 0) = qface[0];
                    QLK(it, jt, kt, 1) = qface[1];
                    QLK(it, jt, kt, 2) = qface[2];
                    QLK(it, jt, kt, 3) = qface[3];
                    QLK(it, jt, kt, 4) = qface[4];
                }
                else
                {
                    //! right boundary
                    QRK(it, jt, kt, 0) = qface[0];
                    QRK(it, jt, kt, 1) = qface[1];
                    QRK(it, jt, kt, 2) = qface[2];
                    QRK(it, jt, kt, 3) = qface[3];
                    QRK(it, jt, kt, 4) = qface[4];
                }
            }
        }
    }


    void NSSolverStruct::CompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int level = grid->GetLevel();
        InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation)
        {
            return;
        }

        int gridtypeofcurrentZone = grid->Type();
        int gridtypeofneighborZone = PHMPI::GetZoneGridType(neighborZoneIndex);
        if (gridtypeofcurrentZone == gridtypeofneighborZone)
        {
            RDouble uploadData = 0.0;
            dataContainer->moveto_begin();
            PHWrite(dataContainer, uploadData);
            return;
        }

        StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

        int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
        int *interfaceIndexContainerForSend = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

        dataContainer->moveto_begin();

        Param_NSSolverStruct *parameters = GetControlParameters();
        int viscousType = parameters->GetViscousType();

        if (viscousType <= INVISCID)
        {
            //! Euler
            RDouble4D &faceInviscidfluxI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxI"));
            RDouble4D &faceInviscidfluxJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxJ"));
            RDouble4D &faceInviscidfluxK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxK"));

            RDouble flux[5];
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForSend[iLocalFace];

                int is, js, ks;
                finestGrid->GetSourceIndexIJK(iFace, 1, is, js, ks);
                finestGrid->RemapMultigridIJK(level, is, js, ks);

                int id, jd, kd;
                int i_lr, j_lr, k_lr;
                finestGrid->GetSourceIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);

                if (i_lr != 0)
                {
                    //! I_face
                    if (i_lr == -1)
                    {
                        //! left boundary
                        flux[0] = faceInviscidfluxI(is, js, ks, 0);
                        flux[1] = faceInviscidfluxI(is, js, ks, 1);
                        flux[2] = faceInviscidfluxI(is, js, ks, 2);
                        flux[3] = faceInviscidfluxI(is, js, ks, 3);
                        flux[4] = faceInviscidfluxI(is, js, ks, 4);
                    }
                    else
                    {
                        //! right boundary
                        flux[0] = -faceInviscidfluxI(is + 1, js, ks, 0);
                        flux[1] = -faceInviscidfluxI(is + 1, js, ks, 1);
                        flux[2] = -faceInviscidfluxI(is + 1, js, ks, 2);
                        flux[3] = -faceInviscidfluxI(is + 1, js, ks, 3);
                        flux[4] = -faceInviscidfluxI(is + 1, js, ks, 4);
                    }
                }
                else if (j_lr != 0)
                {
                    //! J_face
                    if (j_lr == -1)
                    {
                        //! left boundary
                        flux[0] = faceInviscidfluxJ(is, js, ks, 0);
                        flux[1] = faceInviscidfluxJ(is, js, ks, 1);
                        flux[2] = faceInviscidfluxJ(is, js, ks, 2);
                        flux[3] = faceInviscidfluxJ(is, js, ks, 3);
                        flux[4] = faceInviscidfluxJ(is, js, ks, 4);
                    }
                    else
                    {
                        //! right boundary
                        flux[0] = -faceInviscidfluxJ(is, js + 1, ks, 0);
                        flux[1] = -faceInviscidfluxJ(is, js + 1, ks, 1);
                        flux[2] = -faceInviscidfluxJ(is, js + 1, ks, 2);
                        flux[3] = -faceInviscidfluxJ(is, js + 1, ks, 3);
                        flux[4] = -faceInviscidfluxJ(is, js + 1, ks, 4);
                    }
                }
                else
                {
                    //! K_face
                    if (k_lr == -1)
                    {
                        //! left boundary
                        flux[0] = faceInviscidfluxK(is, js, ks, 0);
                        flux[1] = faceInviscidfluxK(is, js, ks, 1);
                        flux[2] = faceInviscidfluxK(is, js, ks, 2);
                        flux[3] = faceInviscidfluxK(is, js, ks, 3);
                        flux[4] = faceInviscidfluxK(is, js, ks, 4);
                    }
                    else
                    {
                        //! right boundary
                        flux[0] = -faceInviscidfluxK(is, js, ks + 1, 0);
                        flux[1] = -faceInviscidfluxK(is, js, ks + 1, 1);
                        flux[2] = -faceInviscidfluxK(is, js, ks + 1, 2);
                        flux[3] = -faceInviscidfluxK(is, js, ks + 1, 3);
                        flux[4] = -faceInviscidfluxK(is, js, ks + 1, 4);
                    }
                }

                for (int m = 0; m < 5; ++ m)
                {
                    PHWrite(dataContainer, flux[m]);
                }
            }
        }
        else
        {
            //! laminar
            RDouble4D &faceInviscidfluxI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxI"));
            RDouble4D &faceInviscidfluxJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxJ"));
            RDouble4D &faceInviscidfluxK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxK"));

            RDouble4D &faceViscousfluxI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxI"));
            RDouble4D &faceViscousfluxJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxJ"));
            RDouble4D &faceViscousfluxK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxK"));

            RDouble Inviscidflux[5];
            RDouble Viscousflux[5];
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForSend[iLocalFace];

                int is, js, ks;
                finestGrid->GetSourceIndexIJK(iFace, 1, is, js, ks);
                finestGrid->RemapMultigridIJK(level, is, js, ks);

                int id, jd, kd;
                int i_lr, j_lr, k_lr;
                finestGrid->GetSourceIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);

                if (i_lr != 0)
                {
                    //! I_face
                    if (i_lr == -1)
                    {
                        //! left boundary
                        Inviscidflux[0] = faceInviscidfluxI(is, js, ks, 0);
                        Inviscidflux[1] = faceInviscidfluxI(is, js, ks, 1);
                        Inviscidflux[2] = faceInviscidfluxI(is, js, ks, 2);
                        Inviscidflux[3] = faceInviscidfluxI(is, js, ks, 3);
                        Inviscidflux[4] = faceInviscidfluxI(is, js, ks, 4);

                        Viscousflux[0] = faceViscousfluxI(is, js, ks, 0);
                        Viscousflux[1] = faceViscousfluxI(is, js, ks, 1);
                        Viscousflux[2] = faceViscousfluxI(is, js, ks, 2);
                        Viscousflux[3] = faceViscousfluxI(is, js, ks, 3);
                        Viscousflux[4] = faceViscousfluxI(is, js, ks, 4);
                    }
                    else
                    {
                        //! right boundary
                        Inviscidflux[0] = -faceInviscidfluxI(is + 1, js, ks, 0);
                        Inviscidflux[1] = -faceInviscidfluxI(is + 1, js, ks, 1);
                        Inviscidflux[2] = -faceInviscidfluxI(is + 1, js, ks, 2);
                        Inviscidflux[3] = -faceInviscidfluxI(is + 1, js, ks, 3);
                        Inviscidflux[4] = -faceInviscidfluxI(is + 1, js, ks, 4);

                        Viscousflux[0] = -faceViscousfluxI(is + 1, js, ks, 0);
                        Viscousflux[1] = -faceViscousfluxI(is + 1, js, ks, 1);
                        Viscousflux[2] = -faceViscousfluxI(is + 1, js, ks, 2);
                        Viscousflux[3] = -faceViscousfluxI(is + 1, js, ks, 3);
                        Viscousflux[4] = -faceViscousfluxI(is + 1, js, ks, 4);
                    }
                }
                else if (j_lr != 0)
                {
                    //! J_face
                    if (j_lr == -1)
                    {
                        //! left boundary
                        Inviscidflux[0] = faceInviscidfluxJ(is, js, ks, 0);
                        Inviscidflux[1] = faceInviscidfluxJ(is, js, ks, 1);
                        Inviscidflux[2] = faceInviscidfluxJ(is, js, ks, 2);
                        Inviscidflux[3] = faceInviscidfluxJ(is, js, ks, 3);
                        Inviscidflux[4] = faceInviscidfluxJ(is, js, ks, 4);

                        Viscousflux[0] = faceViscousfluxJ(is, js, ks, 0);
                        Viscousflux[1] = faceViscousfluxJ(is, js, ks, 1);
                        Viscousflux[2] = faceViscousfluxJ(is, js, ks, 2);
                        Viscousflux[3] = faceViscousfluxJ(is, js, ks, 3);
                        Viscousflux[4] = faceViscousfluxJ(is, js, ks, 4);
                    }
                    else
                    {
                        //! right boundary
                        Inviscidflux[0] = -faceInviscidfluxJ(is, js + 1, ks, 0);
                        Inviscidflux[1] = -faceInviscidfluxJ(is, js + 1, ks, 1);
                        Inviscidflux[2] = -faceInviscidfluxJ(is, js + 1, ks, 2);
                        Inviscidflux[3] = -faceInviscidfluxJ(is, js + 1, ks, 3);
                        Inviscidflux[4] = -faceInviscidfluxJ(is, js + 1, ks, 4);

                        Viscousflux[0] = -faceViscousfluxJ(is, js + 1, ks, 0);
                        Viscousflux[1] = -faceViscousfluxJ(is, js + 1, ks, 1);
                        Viscousflux[2] = -faceViscousfluxJ(is, js + 1, ks, 2);
                        Viscousflux[3] = -faceViscousfluxJ(is, js + 1, ks, 3);
                        Viscousflux[4] = -faceViscousfluxJ(is, js + 1, ks, 4);
                    }
                }
                else
                {
                    //! K_face
                    if (k_lr == -1)
                    {
                        //! left boundary
                        Inviscidflux[0] = faceInviscidfluxK(is, js, ks, 0);
                        Inviscidflux[1] = faceInviscidfluxK(is, js, ks, 1);
                        Inviscidflux[2] = faceInviscidfluxK(is, js, ks, 2);
                        Inviscidflux[3] = faceInviscidfluxK(is, js, ks, 3);
                        Inviscidflux[4] = faceInviscidfluxK(is, js, ks, 4);

                        Viscousflux[0] = faceViscousfluxK(is, js, ks, 0);
                        Viscousflux[1] = faceViscousfluxK(is, js, ks, 1);
                        Viscousflux[2] = faceViscousfluxK(is, js, ks, 2);
                        Viscousflux[3] = faceViscousfluxK(is, js, ks, 3);
                        Viscousflux[4] = faceViscousfluxK(is, js, ks, 4);
                    }
                    else
                    {
                        //! right boundary
                        Inviscidflux[0] = -faceInviscidfluxK(is, js, ks + 1, 0);
                        Inviscidflux[1] = -faceInviscidfluxK(is, js, ks + 1, 1);
                        Inviscidflux[2] = -faceInviscidfluxK(is, js, ks + 1, 2);
                        Inviscidflux[3] = -faceInviscidfluxK(is, js, ks + 1, 3);
                        Inviscidflux[4] = -faceInviscidfluxK(is, js, ks + 1, 4);

                        Viscousflux[0] = -faceViscousfluxK(is, js, ks + 1, 0);
                        Viscousflux[1] = -faceViscousfluxK(is, js, ks + 1, 1);
                        Viscousflux[2] = -faceViscousfluxK(is, js, ks + 1, 2);
                        Viscousflux[3] = -faceViscousfluxK(is, js, ks + 1, 3);
                        Viscousflux[4] = -faceViscousfluxK(is, js, ks + 1, 4);
                    }
                }
                for (int m = 0; m < 5; ++ m)
                {
                    PHWrite(dataContainer, Inviscidflux[m]);
                    PHWrite(dataContainer, Viscousflux[m]);
                }
            }
        }

    }

    void NSSolverStruct::DecompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int level = grid->GetLevel();
        InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation)
        {
            return;
        }

        int gridtypeofcurrentZone = grid->Type();
        int gridtypeofneighborZone = PHMPI::GetZoneGridType(neighborZoneIndex);
        if (gridtypeofcurrentZone == gridtypeofneighborZone)
        {
            RDouble downdData = 0.0;
            dataContainer->moveto_begin();
            PHRead(dataContainer, downdData);
            return;
        }

        StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

        int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
        int *interfaceIndexContainerForReceive = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

        dataContainer->moveto_begin();

        Param_NSSolverStruct *parameters = GetControlParameters();
        int viscousType = parameters->GetViscousType();

        RDouble coefofstrflux = GlobalDataBase::GetDoubleParaFromDB("coefofstrflux");
        RDouble coefofunstrflux = 1.0 - coefofstrflux;

        if (viscousType <= INVISCID)
        {
            //! Euler
            RDouble4D &faceInviscidfluxI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxI"));
            RDouble4D &faceInviscidfluxJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxJ"));
            RDouble4D &faceInviscidfluxK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxK"));

            RDouble flux[5];
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int it, jt, kt;
                finestGrid->GetTargetIndexIJK(iFace, 1, it, jt, kt);
                finestGrid->RemapMultigridIJK(level, it, jt, kt);

                for (int m = 0; m < 5; ++ m)
                {
                    PHRead(dataContainer, flux[m]);
                }

                int id, jd, kd;
                int i_lr, j_lr, k_lr;
                finestGrid->GetTargetIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);
                if (i_lr != 0)
                {
                    //! I_face
                    if (i_lr == -1)
                    {
                        //! left boundary
                        faceInviscidfluxI(it + 1, jt, kt, 0) = (coefofunstrflux * flux[0] + coefofstrflux * faceInviscidfluxI(it + 1, jt, kt, 0));
                        faceInviscidfluxI(it + 1, jt, kt, 1) = (coefofunstrflux * flux[1] + coefofstrflux * faceInviscidfluxI(it + 1, jt, kt, 1));
                        faceInviscidfluxI(it + 1, jt, kt, 2) = (coefofunstrflux * flux[2] + coefofstrflux * faceInviscidfluxI(it + 1, jt, kt, 2));
                        faceInviscidfluxI(it + 1, jt, kt, 3) = (coefofunstrflux * flux[3] + coefofstrflux * faceInviscidfluxI(it + 1, jt, kt, 3));
                        faceInviscidfluxI(it + 1, jt, kt, 4) = (coefofunstrflux * flux[4] + coefofstrflux * faceInviscidfluxI(it + 1, jt, kt, 4));
                    }
                    else
                    {
                        //! right boundary
                        faceInviscidfluxI(it, jt, kt, 0) = (-coefofunstrflux * flux[0] + coefofstrflux * faceInviscidfluxI(it, jt, kt, 0));
                        faceInviscidfluxI(it, jt, kt, 1) = (-coefofunstrflux * flux[1] + coefofstrflux * faceInviscidfluxI(it, jt, kt, 1));
                        faceInviscidfluxI(it, jt, kt, 2) = (-coefofunstrflux * flux[2] + coefofstrflux * faceInviscidfluxI(it, jt, kt, 2));
                        faceInviscidfluxI(it, jt, kt, 3) = (-coefofunstrflux * flux[3] + coefofstrflux * faceInviscidfluxI(it, jt, kt, 3));
                        faceInviscidfluxI(it, jt, kt, 4) = (-coefofunstrflux * flux[4] + coefofstrflux * faceInviscidfluxI(it, jt, kt, 4));
                    }
                }
                else if (j_lr != 0)
                {
                    //! J_face
                    if (j_lr == -1)
                    {
                        //! left boundary
                        faceInviscidfluxJ(it, jt + 1, kt, 0) = (coefofunstrflux * flux[0] + coefofstrflux * faceInviscidfluxJ(it, jt + 1, kt, 0));
                        faceInviscidfluxJ(it, jt + 1, kt, 1) = (coefofunstrflux * flux[1] + coefofstrflux * faceInviscidfluxJ(it, jt + 1, kt, 1));
                        faceInviscidfluxJ(it, jt + 1, kt, 2) = (coefofunstrflux * flux[2] + coefofstrflux * faceInviscidfluxJ(it, jt + 1, kt, 2));
                        faceInviscidfluxJ(it, jt + 1, kt, 3) = (coefofunstrflux * flux[3] + coefofstrflux * faceInviscidfluxJ(it, jt + 1, kt, 3));
                        faceInviscidfluxJ(it, jt + 1, kt, 4) = (coefofunstrflux * flux[4] + coefofstrflux * faceInviscidfluxJ(it, jt + 1, kt, 4));
                    }
                    else
                    {
                        //! right boundary
                        faceInviscidfluxJ(it, jt, kt, 0) = (-coefofunstrflux * flux[0] + coefofstrflux * faceInviscidfluxJ(it, jt, kt, 0));
                        faceInviscidfluxJ(it, jt, kt, 1) = (-coefofunstrflux * flux[1] + coefofstrflux * faceInviscidfluxJ(it, jt, kt, 1));
                        faceInviscidfluxJ(it, jt, kt, 2) = (-coefofunstrflux * flux[2] + coefofstrflux * faceInviscidfluxJ(it, jt, kt, 2));
                        faceInviscidfluxJ(it, jt, kt, 3) = (-coefofunstrflux * flux[3] + coefofstrflux * faceInviscidfluxJ(it, jt, kt, 3));
                        faceInviscidfluxJ(it, jt, kt, 4) = (-coefofunstrflux * flux[4] + coefofstrflux * faceInviscidfluxJ(it, jt, kt, 4));
                    }
                }
                else
                {
                    //! K_face
                    if (k_lr == -1)
                    {
                        //! left boundary
                        faceInviscidfluxK(it, jt, kt + 1, 0) = (coefofunstrflux * flux[0] + coefofstrflux * faceInviscidfluxK(it, jt, kt + 1, 0));
                        faceInviscidfluxK(it, jt, kt + 1, 1) = (coefofunstrflux * flux[1] + coefofstrflux * faceInviscidfluxK(it, jt, kt + 1, 1));
                        faceInviscidfluxK(it, jt, kt + 1, 2) = (coefofunstrflux * flux[2] + coefofstrflux * faceInviscidfluxK(it, jt, kt + 1, 2));
                        faceInviscidfluxK(it, jt, kt + 1, 3) = (coefofunstrflux * flux[3] + coefofstrflux * faceInviscidfluxK(it, jt, kt + 1, 3));
                        faceInviscidfluxK(it, jt, kt + 1, 4) = (coefofunstrflux * flux[4] + coefofstrflux * faceInviscidfluxK(it, jt, kt + 1, 4));
                    }
                    else
                    {
                        //! right boundary
                        faceInviscidfluxK(it, jt, kt, 0) = (-coefofunstrflux * flux[0] + coefofstrflux * faceInviscidfluxK(it, jt, kt, 0));
                        faceInviscidfluxK(it, jt, kt, 1) = (-coefofunstrflux * flux[1] + coefofstrflux * faceInviscidfluxK(it, jt, kt, 1));
                        faceInviscidfluxK(it, jt, kt, 2) = (-coefofunstrflux * flux[2] + coefofstrflux * faceInviscidfluxK(it, jt, kt, 2));
                        faceInviscidfluxK(it, jt, kt, 3) = (-coefofunstrflux * flux[3] + coefofstrflux * faceInviscidfluxK(it, jt, kt, 3));
                        faceInviscidfluxK(it, jt, kt, 4) = (-coefofunstrflux * flux[4] + coefofstrflux * faceInviscidfluxK(it, jt, kt, 4));
                    }
                }
            }
        }
        else
        {
            //! laminar
            RDouble4D &faceInviscidfluxI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxI"));
            RDouble4D &faceInviscidfluxJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxJ"));
            RDouble4D &faceInviscidfluxK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceInviscidfluxK"));

            RDouble4D &faceViscousfluxI = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxI"));
            RDouble4D &faceViscousfluxJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxJ"));
            RDouble4D &faceViscousfluxK = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("faceViscousfluxK"));

            RDouble Inviscidflux[5];
            RDouble Viscousflux[5];
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int it, jt, kt;
                finestGrid->GetTargetIndexIJK(iFace, 1, it, jt, kt);
                finestGrid->RemapMultigridIJK(level, it, jt, kt);

                for (int m = 0; m < 5; ++ m)
                {
                    PHRead(dataContainer, Inviscidflux[m]);
                    PHRead(dataContainer, Viscousflux[m]);
                }

                int id, jd, kd;
                int i_lr, j_lr, k_lr;
                finestGrid->GetTargetIndexIJK_fornodetransfer(iFace, id, jd, kd, i_lr, j_lr, k_lr);
                if (i_lr != 0)
                {
                    //! I_face
                    if (i_lr == -1)
                    {
                        //! left boundary
                        faceInviscidfluxI(it + 1, jt, kt, 0) = (coefofunstrflux * Inviscidflux[0] + coefofstrflux * faceInviscidfluxI(it + 1, jt, kt, 0));
                        faceInviscidfluxI(it + 1, jt, kt, 1) = (coefofunstrflux * Inviscidflux[1] + coefofstrflux * faceInviscidfluxI(it + 1, jt, kt, 1));
                        faceInviscidfluxI(it + 1, jt, kt, 2) = (coefofunstrflux * Inviscidflux[2] + coefofstrflux * faceInviscidfluxI(it + 1, jt, kt, 2));
                        faceInviscidfluxI(it + 1, jt, kt, 3) = (coefofunstrflux * Inviscidflux[3] + coefofstrflux * faceInviscidfluxI(it + 1, jt, kt, 3));
                        faceInviscidfluxI(it + 1, jt, kt, 4) = (coefofunstrflux * Inviscidflux[4] + coefofstrflux * faceInviscidfluxI(it + 1, jt, kt, 4));

                        faceViscousfluxI(it + 1, jt, kt, 0) = (coefofunstrflux * Viscousflux[0] + coefofstrflux * faceViscousfluxI(it + 1, jt, kt, 0));
                        faceViscousfluxI(it + 1, jt, kt, 1) = (coefofunstrflux * Viscousflux[1] + coefofstrflux * faceViscousfluxI(it + 1, jt, kt, 1));
                        faceViscousfluxI(it + 1, jt, kt, 2) = (coefofunstrflux * Viscousflux[2] + coefofstrflux * faceViscousfluxI(it + 1, jt, kt, 2));
                        faceViscousfluxI(it + 1, jt, kt, 3) = (coefofunstrflux * Viscousflux[3] + coefofstrflux * faceViscousfluxI(it + 1, jt, kt, 3));
                        faceViscousfluxI(it + 1, jt, kt, 4) = (coefofunstrflux * Viscousflux[4] + coefofstrflux * faceViscousfluxI(it + 1, jt, kt, 4));
                    }
                    else
                    {
                        //! right boundary
                        faceInviscidfluxI(it, jt, kt, 0) = (-coefofunstrflux * Inviscidflux[0] + coefofstrflux * faceInviscidfluxI(it, jt, kt, 0));
                        faceInviscidfluxI(it, jt, kt, 1) = (-coefofunstrflux * Inviscidflux[1] + coefofstrflux * faceInviscidfluxI(it, jt, kt, 1));
                        faceInviscidfluxI(it, jt, kt, 2) = (-coefofunstrflux * Inviscidflux[2] + coefofstrflux * faceInviscidfluxI(it, jt, kt, 2));
                        faceInviscidfluxI(it, jt, kt, 3) = (-coefofunstrflux * Inviscidflux[3] + coefofstrflux * faceInviscidfluxI(it, jt, kt, 3));
                        faceInviscidfluxI(it, jt, kt, 4) = (-coefofunstrflux * Inviscidflux[4] + coefofstrflux * faceInviscidfluxI(it, jt, kt, 4));

                        faceViscousfluxI(it, jt, kt, 0) = (-coefofunstrflux * Viscousflux[0] + coefofstrflux * faceViscousfluxI(it, jt, kt, 0));
                        faceViscousfluxI(it, jt, kt, 1) = (-coefofunstrflux * Viscousflux[1] + coefofstrflux * faceViscousfluxI(it, jt, kt, 1));
                        faceViscousfluxI(it, jt, kt, 2) = (-coefofunstrflux * Viscousflux[2] + coefofstrflux * faceViscousfluxI(it, jt, kt, 2));
                        faceViscousfluxI(it, jt, kt, 3) = (-coefofunstrflux * Viscousflux[3] + coefofstrflux * faceViscousfluxI(it, jt, kt, 3));
                        faceViscousfluxI(it, jt, kt, 4) = (-coefofunstrflux * Viscousflux[4] + coefofstrflux * faceViscousfluxI(it, jt, kt, 4));
                    }
                }
                else if (j_lr != 0)
                {
                    //! J_face
                    if (j_lr == -1)
                    {
                        //! left boundary
                        faceInviscidfluxJ(it, jt + 1, kt, 0) = (coefofunstrflux * Inviscidflux[0] + coefofstrflux * faceInviscidfluxJ(it, jt + 1, kt, 0));
                        faceInviscidfluxJ(it, jt + 1, kt, 1) = (coefofunstrflux * Inviscidflux[1] + coefofstrflux * faceInviscidfluxJ(it, jt + 1, kt, 1));
                        faceInviscidfluxJ(it, jt + 1, kt, 2) = (coefofunstrflux * Inviscidflux[2] + coefofstrflux * faceInviscidfluxJ(it, jt + 1, kt, 2));
                        faceInviscidfluxJ(it, jt + 1, kt, 3) = (coefofunstrflux * Inviscidflux[3] + coefofstrflux * faceInviscidfluxJ(it, jt + 1, kt, 3));
                        faceInviscidfluxJ(it, jt + 1, kt, 4) = (coefofunstrflux * Inviscidflux[4] + coefofstrflux * faceInviscidfluxJ(it, jt + 1, kt, 4));

                        faceViscousfluxJ(it, jt + 1, kt, 0) = (coefofunstrflux * Viscousflux[0] + coefofstrflux * faceViscousfluxJ(it, jt + 1, kt, 0));
                        faceViscousfluxJ(it, jt + 1, kt, 1) = (coefofunstrflux * Viscousflux[1] + coefofstrflux * faceViscousfluxJ(it, jt + 1, kt, 1));
                        faceViscousfluxJ(it, jt + 1, kt, 2) = (coefofunstrflux * Viscousflux[2] + coefofstrflux * faceViscousfluxJ(it, jt + 1, kt, 2));
                        faceViscousfluxJ(it, jt + 1, kt, 3) = (coefofunstrflux * Viscousflux[3] + coefofstrflux * faceViscousfluxJ(it, jt + 1, kt, 3));
                        faceViscousfluxJ(it, jt + 1, kt, 4) = (coefofunstrflux * Viscousflux[4] + coefofstrflux * faceViscousfluxJ(it, jt + 1, kt, 4));
                    }
                    else
                    {
                        //! right boundary
                        faceInviscidfluxJ(it, jt, kt, 0) = (-coefofunstrflux * Inviscidflux[0] + coefofstrflux * faceInviscidfluxJ(it, jt, kt, 0));
                        faceInviscidfluxJ(it, jt, kt, 1) = (-coefofunstrflux * Inviscidflux[1] + coefofstrflux * faceInviscidfluxJ(it, jt, kt, 1));
                        faceInviscidfluxJ(it, jt, kt, 2) = (-coefofunstrflux * Inviscidflux[2] + coefofstrflux * faceInviscidfluxJ(it, jt, kt, 2));
                        faceInviscidfluxJ(it, jt, kt, 3) = (-coefofunstrflux * Inviscidflux[3] + coefofstrflux * faceInviscidfluxJ(it, jt, kt, 3));
                        faceInviscidfluxJ(it, jt, kt, 4) = (-coefofunstrflux * Inviscidflux[4] + coefofstrflux * faceInviscidfluxJ(it, jt, kt, 4));

                        faceViscousfluxJ(it, jt, kt, 0) = (-coefofunstrflux * Viscousflux[0] + coefofstrflux * faceViscousfluxJ(it, jt, kt, 0));
                        faceViscousfluxJ(it, jt, kt, 1) = (-coefofunstrflux * Viscousflux[1] + coefofstrflux * faceViscousfluxJ(it, jt, kt, 1));
                        faceViscousfluxJ(it, jt, kt, 2) = (-coefofunstrflux * Viscousflux[2] + coefofstrflux * faceViscousfluxJ(it, jt, kt, 2));
                        faceViscousfluxJ(it, jt, kt, 3) = (-coefofunstrflux * Viscousflux[3] + coefofstrflux * faceViscousfluxJ(it, jt, kt, 3));
                        faceViscousfluxJ(it, jt, kt, 4) = (-coefofunstrflux * Viscousflux[4] + coefofstrflux * faceViscousfluxJ(it, jt, kt, 4));
                    }
                }
                else
                {
                    //! K_face
                    if (k_lr == -1)
                    {
                        //! left boundary
                        faceInviscidfluxK(it, jt, kt + 1, 0) = (coefofunstrflux * Inviscidflux[0] + coefofstrflux * faceInviscidfluxK(it, jt, kt + 1, 0));
                        faceInviscidfluxK(it, jt, kt + 1, 1) = (coefofunstrflux * Inviscidflux[1] + coefofstrflux * faceInviscidfluxK(it, jt, kt + 1, 1));
                        faceInviscidfluxK(it, jt, kt + 1, 2) = (coefofunstrflux * Inviscidflux[2] + coefofstrflux * faceInviscidfluxK(it, jt, kt + 1, 2));
                        faceInviscidfluxK(it, jt, kt + 1, 3) = (coefofunstrflux * Inviscidflux[3] + coefofstrflux * faceInviscidfluxK(it, jt, kt + 1, 3));
                        faceInviscidfluxK(it, jt, kt + 1, 4) = (coefofunstrflux * Inviscidflux[4] + coefofstrflux * faceInviscidfluxK(it, jt, kt + 1, 4));

                        faceViscousfluxK(it, jt, kt + 1, 0) = (coefofunstrflux * Viscousflux[0] + coefofstrflux * faceViscousfluxK(it, jt, kt + 1, 0));
                        faceViscousfluxK(it, jt, kt + 1, 1) = (coefofunstrflux * Viscousflux[1] + coefofstrflux * faceViscousfluxK(it, jt, kt + 1, 1));
                        faceViscousfluxK(it, jt, kt + 1, 2) = (coefofunstrflux * Viscousflux[2] + coefofstrflux * faceViscousfluxK(it, jt, kt + 1, 2));
                        faceViscousfluxK(it, jt, kt + 1, 3) = (coefofunstrflux * Viscousflux[3] + coefofstrflux * faceViscousfluxK(it, jt, kt + 1, 3));
                        faceViscousfluxK(it, jt, kt + 1, 4) = (coefofunstrflux * Viscousflux[4] + coefofstrflux * faceViscousfluxK(it, jt, kt + 1, 4));
                    }
                    else
                    {
                        //! right boundary
                        faceInviscidfluxK(it, jt, kt, 0) = (-coefofunstrflux * Inviscidflux[0] + coefofstrflux * faceInviscidfluxK(it, jt, kt, 0));
                        faceInviscidfluxK(it, jt, kt, 1) = (-coefofunstrflux * Inviscidflux[1] + coefofstrflux * faceInviscidfluxK(it, jt, kt, 1));
                        faceInviscidfluxK(it, jt, kt, 2) = (-coefofunstrflux * Inviscidflux[2] + coefofstrflux * faceInviscidfluxK(it, jt, kt, 2));
                        faceInviscidfluxK(it, jt, kt, 3) = (-coefofunstrflux * Inviscidflux[3] + coefofstrflux * faceInviscidfluxK(it, jt, kt, 3));
                        faceInviscidfluxK(it, jt, kt, 4) = (-coefofunstrflux * Inviscidflux[4] + coefofstrflux * faceInviscidfluxK(it, jt, kt, 4));

                        faceViscousfluxK(it, jt, kt, 0) = (-coefofunstrflux * Viscousflux[0] + coefofstrflux * faceViscousfluxK(it, jt, kt, 0));
                        faceViscousfluxK(it, jt, kt, 1) = (-coefofunstrflux * Viscousflux[1] + coefofstrflux * faceViscousfluxK(it, jt, kt, 1));
                        faceViscousfluxK(it, jt, kt, 2) = (-coefofunstrflux * Viscousflux[2] + coefofstrflux * faceViscousfluxK(it, jt, kt, 2));
                        faceViscousfluxK(it, jt, kt, 3) = (-coefofunstrflux * Viscousflux[3] + coefofstrflux * faceViscousfluxK(it, jt, kt, 3));
                        faceViscousfluxK(it, jt, kt, 4) = (-coefofunstrflux * Viscousflux[4] + coefofstrflux * faceViscousfluxK(it, jt, kt, 4));
                    }
                }
            }
        }
    }

    void NSSolverStruct::CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int level = grid->GetLevel();
        InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation)
        {
            return;
        }

        StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

        int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
        int *interfaceIndexContainerForSend = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

        dataContainer->moveto_begin();

        RDouble4D &fieldSend = fieldProxy->GetField_STR();

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForSend[iLocalFace];
                int is, js, ks;
                finestGrid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
                finestGrid->RemapMultigridIJK(level, is, js, ks);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHWrite(dataContainer, fieldSend(is, js, ks, m));
                }
            }
        }
    }

    void NSSolverStruct::DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int level = grid->GetLevel();
        InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
        if (!interfaceInformation)
        {
            return;
        }

        StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
        InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();
        CompositeStructBC *composite_bcregion = finestGrid->GetCompositeBCRegion();

        int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
        int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
        int *interfaceIndexContainerForReceive = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

        int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
        double rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
        rotationAngle = rotationAngle * PI / 180;

        dataContainer->moveto_begin();

        RDouble4D &fieldRecv = fieldProxy->GetField_STR();

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int it, jt, kt;
                finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
                finestGrid->RemapMultigridIJK(level, it, jt, kt);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHRead(dataContainer, fieldRecv(it, jt, kt, m));
                }
                if (periodicType == 2)
                {
                    int *ibcregions = composite_bcregion->GetIFaceInfo();
                    int iBCRegion = ibcregions[iFace];
                    StructBC *bcregion = composite_bcregion->GetBCRegion(iBCRegion);
                    string bcName = bcregion->GetBCName();
                    if (bcName == "Periodic_up" || bcName == "Periodic_down")
                    {
                        double rotfg[3] = { 0, 0, 0 };
                        if (bcName == "Periodic_up")
                        {
                            rotfg[0] = fieldRecv(it, jt, kt, 1);
                            rotfg[1] = fieldRecv(it, jt, kt, 2) * cos(2 * PI - rotationAngle) - fieldRecv(it, jt, kt, 3) * sin(2 * PI - rotationAngle);
                            rotfg[2] = fieldRecv(it, jt, kt, 2) * sin(2 * PI - rotationAngle) + fieldRecv(it, jt, kt, 3) * cos(2 * PI - rotationAngle);

                        }
                        else if (bcName == "Periodic_down")
                        {
                            rotfg[0] = fieldRecv(it, jt, kt, 1);
                            rotfg[1] = fieldRecv(it, jt, kt, 2) * cos(rotationAngle) - fieldRecv(it, jt, kt, 3) * sin(rotationAngle);
                            rotfg[2] = fieldRecv(it, jt, kt, 2) * sin(rotationAngle) + fieldRecv(it, jt, kt, 3) * cos(rotationAngle);
                        }

                        for (int m = 1; m <= 3; ++ m)
                        {
                            fieldRecv(it, jt, kt, m) = rotfg[m - 1];
                        }
                    }
                }
            }
        }
    }

    GradientCellCenter *NSSolverStruct::GetGradientCellCenterUVWT(Grid *gridIn, FieldProxy *primitiveVariables)
    {
        int level = gridIn->GetLevel();
        if (gradientCellCenterUVWT[level] == NULL)
        {
            gradientCellCenterUVWT[level] = new GradientCellCenter(gridIn, STRUCTGRID, primitiveVariables, 4);
        }
        return gradientCellCenterUVWT[level];
    }

    GradientCellCenter *NSSolverStruct::GetGradientCellCenterUVWT(Grid *gridIn)
    {
        int level = gridIn->GetLevel();
        return gradientCellCenterUVWT[level];
    }

    LIB_EXPORT void NSSolverStruct::ComputePostVisualVariables(Post_Visual *postVisualization)
    {
        StructGrid *grid = StructGridCast(postVisualization->GetGrid());


        if (postVisualization->GetFlowType() == AverageFlow)
        {
            RDouble4D &qAverage = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage"));

            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DENSITY);
            postVisualization->UpdateVisualNodeVarPtr(varName, qAverage, IDX::IR);

            varName = postVisualization->GetVariableName(PHSPACE::VISUAL_U);
            postVisualization->UpdateVisualNodeVarPtr(varName, qAverage, IDX::IU);

            varName = postVisualization->GetVariableName(PHSPACE::VISUAL_V);
            postVisualization->UpdateVisualNodeVarPtr(varName, qAverage, IDX::IV);

            varName = postVisualization->GetVariableName(PHSPACE::VISUAL_W);
            postVisualization->UpdateVisualNodeVarPtr(varName, qAverage, IDX::IW);

            varName = postVisualization->GetVariableName(PHSPACE::VISUAL_PRESSURE);
            postVisualization->UpdateVisualNodeVarPtr(varName, qAverage, IDX::IP);

            return;
        }

        //! added by zzp 202108, for ReynoldsStress output
        if (postVisualization->GetFlowType() == AverageReynoldsStress)
        {
            RDouble4D &tauAverage = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("tauAverage"));

            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_XX);
            postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage, 0);

            varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_YY);
            postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage, 1);

            varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_ZZ);
            postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage, 2);

            varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_XY);
            postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage, 3);

            varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_XZ);
            postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage, 4);

            varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_YZ);
            postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage, 5);

            return;
        }

        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        int isOverset = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
        if (isOverset)
        {
            int ni = grid->GetNI();
            int nj = grid->GetNJ();
            int nk = grid->GetNK();

            Range II(1, ni - 1);
            Range JJ(1, nj - 1);
            Range KK(1, nk - 1);
            if (nk == 1) KK.setRange(1, 1);

            Int3D &iblank = *new Int3D(II, JJ, KK, fortranArray);
            int *iBlank = UnstructGridCast(grid)->GetBlankIndex();
            int numberOfCells = grid->GetNTotalCell();
            for (int iCell = 0; iCell < numberOfCells; ++ iCell)
            {
                int i, j, k;
                DecodeIJK(iCell, i, j, k, ni - 1, nj - 1, nk - 1);
                iblank(i, j, k) = iBlank[iCell];
            }

            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_IBLANK))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_IBLANK);
                postVisualization->UpdateVisualNodeVarPtr(varName, iblank);
            }
            delete &iblank;
        }
        Param_NSSolverStruct *parameters = GetControlParameters();

        using namespace IDX;
        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DENSITY))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DENSITY);
            postVisualization->UpdateVisualNodeVarPtr(varName, primitiveVariables, IDX::IR);
        }
        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_U))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_U);
            postVisualization->UpdateVisualNodeVarPtr(varName, primitiveVariables, IDX::IU);
        }
        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_V))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_V);
            postVisualization->UpdateVisualNodeVarPtr(varName, primitiveVariables, IDX::IV);
        }
        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_W))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_W);
            postVisualization->UpdateVisualNodeVarPtr(varName, primitiveVariables, IDX::IW);
        }
        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_PRESSURE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_PRESSURE);
            postVisualization->UpdateVisualNodeVarPtr(varName, primitiveVariables, IDX::IP);
        }

        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TEMPERATURE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TEMPERATURE);
            postVisualization->UpdateVisualNodeVarPtr(varName, temperatures, IDX::ITT);
        }

        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MACH))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MACH);

            int isOverset = parameters->GetIsOverLapping();

            RDouble4D *mach = 0;
            if (!isOverset)
            {
                mach = CompMachNumber(grid);
            }
            else
            {
                mach = ComputeMachNumberField(grid);
            }

            postVisualization->UpdateVisualNodeVarPtr(varName, *mach, 0);
            delete mach;
        }

        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STREAMLINE_U))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STREAMLINE_U);
            postVisualization->UpdateVisualNodeVarPtr(varName, primitiveVariables, IDX::IU);
        }

        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STREAMLINE_V))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STREAMLINE_V);
            postVisualization->UpdateVisualNodeVarPtr(varName, primitiveVariables, IDX::IV);
        }

        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STREAMLINE_W))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STREAMLINE_W);
            postVisualization->UpdateVisualNodeVarPtr(varName, primitiveVariables, IDX::IW);
        }

        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TIME_STEP))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TIME_STEP);
            RDouble3D &dt = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("dt"));
            postVisualization->UpdateVisualNodeVarPtr(varName, dt);
        }
        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VOLUME))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VOLUME);
            RDouble3D &vol = *(grid->GetCellVolume());
            postVisualization->UpdateVisualNodeVarPtr(varName, vol);
        }
        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_WALL_DIST))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_WALL_DIST);
            RDouble3D &wallDist = *(grid->GetWallDist());
            postVisualization->UpdateVisualNodeVarPtr(varName, wallDist);
        }

        //! Bell 20130319 add.
        int viscousType = parameters->GetViscousType();
        string viscousName = parameters->GetViscousName();
        int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
        int transitionType = GlobalDataBase::GetIntParaFromDB("transitionType");

        if (viscousType == LAMINAR && iLES == NOLES_SOLVER)
        {
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_LAMINAR))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_LAMINAR);
                RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
                postVisualization->UpdateVisualNodeVarPtr(varName, viscousLaminar);
            }
        }
        else if (viscousType > LAMINAR || iLES == LES_SOLVER)
        {
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_LAMINAR))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_LAMINAR);
                RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
                postVisualization->UpdateVisualNodeVarPtr(varName, viscousLaminar);
            }

            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_TURBULENT))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_TURBULENT);
                RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
                postVisualization->UpdateVisualNodeVarPtr(varName, viscousTurbulence);
            }
            //! For two-equation output turbulent kinetic energy and dissipative rate.
            if (viscousName.substr(0, 6) == "2eq-kw")
            {
                RDouble4D &q_turb = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q_turb"));
                if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MODELED_TKE))
                {
                    string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MODELED_TKE);
                    postVisualization->UpdateVisualNodeVarPtr(varName, q_turb, IDX::IKE);
                }
                if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MODELED_DISSIPATION))
                {
                    string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MODELED_DISSIPATION);
                    postVisualization->UpdateVisualNodeVarPtr(varName, q_turb, IDX::IKW);
                }
                if (viscousName.substr(0, 17) == "2eq-kw-menter-sst")   //! for blending function F1 and F2.
                {
                    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_SST_F1))
                    {
                        RDouble3D &blend = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("blend"));   // for blengding function F1
                        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_SST_F1);
                        postVisualization->UpdateVisualNodeVarPtr(varName, blend);
                    }
                    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_SST_F2))
                    {
                        RDouble3D &SST_F2 = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("SST_F2"));   //! for blengding function F2
                        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_SST_F2);
                        postVisualization->UpdateVisualNodeVarPtr(varName, SST_F2);
                    }
                    if (transitionType == IREGAMA)
                    {
                        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_INTERMITTENCY))
                        {
                            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_INTERMITTENCY);
                            postVisualization->UpdateVisualNodeVarPtr(varName, q_turb, IDX::IGAMA);
                        }
                        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MOMENTUMTHICKREYNOLDS))
                        {
                            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MOMENTUMTHICKREYNOLDS);
                            postVisualization->UpdateVisualNodeVarPtr(varName, q_turb, IDX::IRECT);
                        }
                        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TRANSITION_GAMAEFF))
                        {
                            RDouble3D &gamaeff = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gamaeff"));   // for transition gamaeff
                            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TRANSITION_GAMAEFF);
                            postVisualization->UpdateVisualNodeVarPtr(varName, gamaeff);
                        }
                        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TRANSITION_RESCF))
                        {
                            RDouble3D &rescf = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("rescf"));   // for transition Rescf
                            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TRANSITION_RESCF);
                            postVisualization->UpdateVisualNodeVarPtr(varName, rescf);
                        }
                    }
                }
            }
        }

        int visualField = parameters->GetPlotFieldType();
        if (WantVisualField(grid) || visualField == TEC_SPACE::BlockVisual)
        {
            //! Vorticity, Vorticity, by Q value.
            RDouble3D *vorticity_x, *vorticity_y, *vorticity_z, *vorticityMagnitude, *strain_rate, *Q_criteria;

            int ni = grid->GetNI();
            int nj = grid->GetNJ();
            int nk = grid->GetNK();

            Range I, J, K;
            GetRange(ni, nj, nk, -2, 1, I, J, K);

            vorticity_x = new RDouble3D(I, J, K, fortranArray);
            vorticity_y = new RDouble3D(I, J, K, fortranArray);
            vorticity_z = new RDouble3D(I, J, K, fortranArray);
            vorticityMagnitude = new RDouble3D(I, J, K, fortranArray);
            strain_rate = new RDouble3D(I, J, K, fortranArray);
            Q_criteria = new RDouble3D(I, J, K, fortranArray);

            ComputeVorticitybyQCriteria(grid, vorticity_x, vorticity_y, vorticity_z, vorticityMagnitude, strain_rate, Q_criteria);

            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_X))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_X);
                postVisualization->UpdateVisualNodeVarPtr(varName, *vorticity_x);
            }

            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_Y))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_Y);
                postVisualization->UpdateVisualNodeVarPtr(varName, *vorticity_y);
            }

            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_Z))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_Z);
                postVisualization->UpdateVisualNodeVarPtr(varName, *vorticity_z);
            }

            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_MAGNITUDE))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_MAGNITUDE);
                postVisualization->UpdateVisualNodeVarPtr(varName, *vorticityMagnitude);
            }

            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STRAIN_RATE))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STRAIN_RATE);
                postVisualization->UpdateVisualNodeVarPtr(varName, *strain_rate);
            }

            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_Q_CRITERIA))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_Q_CRITERIA);
                postVisualization->UpdateVisualNodeVarPtr(varName, *Q_criteria);
            }

            delete vorticity_x;
            delete vorticity_y;
            delete vorticity_z;
            delete vorticityMagnitude;
            delete strain_rate;
            delete Q_criteria;

        }
        else if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_X)
            || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_Y)
            || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_Z)
            || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_MAGNITUDE)
            || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STRAIN_RATE)
            || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_Q_CRITERIA))
        {
            TK_Exit::ExceptionExit("Error: can not visualization vorticity when plotFieldType = 0!");
        }

        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_CP))
        {
            RDouble4D *cp;
            ComputeCp(grid, cp);
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_CP);
            postVisualization->UpdateVisualNodeVarPtr(varName, *cp, 0);
            delete cp;
        }

        //! Chemical species must be the last!
        //using namespace GAS_SPACE;
        int nNSEquation = parameters->GetNSEquationNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();
        int numberOfSpecies = parameters->GetNumberOfSpecies();

        //! Write vibration temperature.
        if (nChemical > 0 && nTemperatureModel > 1)
        {
            //To write the vibratial temperature and vibrational energy.
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TEMPERATURE_VIBRATION))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TEMPERATURE_VIBRATION);
                postVisualization->UpdateVisualNodeVarPtr(varName, temperatures, IDX::ITV);
            }
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_ENERGY_VIBRATION))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_ENERGY_VIBRATION);
                postVisualization->UpdateVisualNodeVarPtr(varName, primitiveVariables, nNSEquation + numberOfSpecies);
            }

            //To write the electric temperature and vibrational energy.
            if (nTemperatureModel == 3)
            {
                if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TEMPERATURE_ELECTRON))
                {
                    string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TEMPERATURE_ELECTRON);
                    postVisualization->UpdateVisualNodeVarPtr(varName, temperatures, IDX::ITE);
                }
                if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_ENERGY_ELECTRON))
                {
                    string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_ENERGY_ELECTRON);
                    postVisualization->UpdateVisualNodeVarPtr(varName, primitiveVariables, nNSEquation + numberOfSpecies + 1);
                }
            }
        }

        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_DENSITY)
            || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_U)
            || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_V)
            || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_W)
            || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_PRESSURE)
            || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_TEMPERATURE))
        {
            RDouble4D *dimensionalVariables = ComputeDimensionalVariables(grid);
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_DENSITY))
            {
                //Write the dimensional value of density.
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_DENSITY);
                postVisualization->UpdateVisualNodeVarPtr(varName, *dimensionalVariables, 0);
            }
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_U))
            {
                //Write the dimensional value of U.
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_U);
                postVisualization->UpdateVisualNodeVarPtr(varName, *dimensionalVariables, 1);
            }
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_V))
            {
                //Write the dimensional value of V.
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_V);
                postVisualization->UpdateVisualNodeVarPtr(varName, *dimensionalVariables, 2);
            }
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_W))
            {
                //Write the dimensional value of W.
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_W);
                postVisualization->UpdateVisualNodeVarPtr(varName, *dimensionalVariables, 3);
            }
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_PRESSURE))
            {
                //Write the dimensional value of pressure.
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_PRESSURE);
                postVisualization->UpdateVisualNodeVarPtr(varName, *dimensionalVariables, 4);
            }
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_TEMPERATURE))
            {
                //Write the dimensional value of temperature.
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_TEMPERATURE);
                postVisualization->UpdateVisualNodeVarPtr(varName, *dimensionalVariables, 5);
            }
            delete dimensionalVariables;
        }

        //! To write the mass fractions of species.
        if (nChemical > 0)
        {
            //Write the mass fractions of species.
            using namespace GAS_SPACE;
            string *varname = gas->GetNameOfSpecies();
            for (int m = 0; m < numberOfSpecies; ++ m)
            {
                string varName = varname[m];
                postVisualization->UpdateVisualNodeVarPtr(varName, primitiveVariables, nNSEquation + m);
            }

            //To write the number density of electron.
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_ELECTRON_NUMBER))
            {
                RDouble4D *tmpValue = ComputeDimensionalElectron(grid);
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_ELECTRON_NUMBER);
                postVisualization->UpdateVisualNodeVarPtr(varName, *tmpValue, 0);
                delete tmpValue;
            }
        }

        if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GAMA))
        {
            int ni = grid->GetNI();
            int nj = grid->GetNJ();
            int nk = grid->GetNK();

            Range I, J, K;
            GetRange(ni, nj, nk, -2, 1, I, J, K);
            Range M(0, 0);
            RDouble4D *refgama = new RDouble4D(I, J, K, M, fortranArray);

            int ist, ied, jst, jed, kst, ked;
            grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

            int nEquation = GetNumberOfEquations();
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_GAMA);
            RDouble3D &gama = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("gama"));
            RDouble *primitiveVars = new RDouble[nEquation];

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            primitiveVars[m] = primitiveVariables(i, j, k, m);
                        }
                        gas->GetSpecificHeatRatio(primitiveVars, gama(i, j, k));
                        (*refgama)(i, j, k, 0) = gama(i, j, k);
                    }
                }
            }
            postVisualization->UpdateVisualNodeVarPtr(varName, *refgama, 0);
            delete refgama;
            delete[] primitiveVars;
        }
    }

    LIB_EXPORT void NSSolverStruct::ComputePostProbesVariables(Post_Probes *postProbesVar)
    {
        StructGrid *grid = StructGridCast(postProbesVar->GetGrid());

        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperatures = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));

        Param_NSSolverStruct *parameters = GetControlParameters();

        using namespace IDX;
        if (postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DENSITY))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_DENSITY);
            postProbesVar->UpdateProbesVarPtr(varName, primitiveVariables, IDX::IR);
        }

        if (postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_U))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_U);
            postProbesVar->UpdateProbesVarPtr(varName, primitiveVariables, IDX::IU);
        }

        if (postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_V))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_V);
            postProbesVar->UpdateProbesVarPtr(varName, primitiveVariables, IDX::IV);
        }

        if (postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_W))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_W);
            postProbesVar->UpdateProbesVarPtr(varName, primitiveVariables, IDX::IW);
        }

        if (postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_PRESSURE))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_PRESSURE);
            postProbesVar->UpdateProbesVarPtr(varName, primitiveVariables, IDX::IP);
        }

        if (postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_TEMPERATURE))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_TEMPERATURE);
            postProbesVar->UpdateProbesVarPtr(varName, temperatures, IDX::ITT);
        }

        if (postProbesVar->IsNeedMonitoring(PHSPACE::VISUAL_MACH))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::VISUAL_MACH);

            int isOverset = parameters->GetIsOverLapping();

            RDouble4D *mach = 0;
            if (!isOverset)
            {
                mach = CompMachNumber(grid);
            }
            else
            {
                mach = ComputeMachNumberField(grid);
            }

            postProbesVar->UpdateProbesVarPtr(varName, *mach, 0);
            delete mach;
        }
    }

    LIB_EXPORT void NSSolverStruct::InitControlParameters()
    {
        FreeControlParameters();
        controlParameters = new Param_NSSolverStruct();
        controlParameters->Init();
    }

    void NSSolverStruct::Monitor(Grid *gridIn, RDouble4D *dq)
    {
        StructGrid *grid = StructGridCast(gridIn);

        int ist, ied, jst, jed, kst, ked;
        grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

        RDouble4D &residual = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("res"));
        RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        RDouble4D &temperature = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("t"));
        RDouble2D *surfaceTemperature = reinterpret_cast <RDouble2D *> (gridIn->GetDataPtr("surfaceTemperature"));
        RDouble4D &diagonal = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("diagonal"));
        int nStep = GlobalDataBase::GetIntParaFromDB("outnstep");
        int nEqn = GetNumberOfEquations();
        ostringstream oss;
        using namespace IDX;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    if (i == 1 && j == 1 && k == 1)
                    {
                        oss << "I=" << i << ", J=" << j << ", K=" << k << endl;
                        oss << nStep;
                        for (int s = 0; s < nEqn; ++ s)
                        {
                            oss << ", " << primitiveVariables(i, j, k, s);
                            //oss<<", "<<residual(i, j, k, s);
                            //oss<<", "<<(*dq)(i, j, k, s);
                            //oss<<", "<<diagonal(i, j, k, s);
                        }
                        oss << endl;
                        oss << ", " << temperature(i, j, k, ITT) << endl;
                    }
                }
            }
        }
        oss << (*surfaceTemperature)(0, 0) << endl;
        WriteLogFile(oss);
    }

    LIB_EXPORT Param_NSSolverStruct *NSSolverStruct::GetControlParameters() const
    {
        return static_cast <Param_NSSolverStruct *> (controlParameters);
    }

}
