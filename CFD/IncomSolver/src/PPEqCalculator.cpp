#include "PPEqCalculator.h"




namespace PHSPACE
{
PPEqCalculator::PPEqCalculator():IncomCalculator()
{
}

PPEqCalculator::~PPEqCalculator()
{

}

void PPEqCalculator::correctPressureAndVelocity(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble *vol = grid->GetCellVolume();
    RDouble PPEqRelaxCoeff = GlobalDataBase::GetDoubleParaFromDB("PPEqRelaxCoeff");
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
    RDouble *dpcdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpcdx"));
    RDouble *dpcdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpcdy"));
    RDouble *dpcdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpcdz"));
    RDouble *pc = reinterpret_cast<RDouble *>(grid->GetDataPtr("pc"));
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble *huU = reinterpret_cast<RDouble *>(grid->GetDataPtr("huU"));
    RDouble *huV = reinterpret_cast<RDouble *>(grid->GetDataPtr("huV"));
    RDouble *huW = reinterpret_cast<RDouble *>(grid->GetDataPtr("huW"));
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY
            || bcType == PHENGLEI::PRESSURE_INLET || bcType == PHENGLEI::INFLOW)
        {
            int *leftCellOfFace = grid->GetLeftCellOfFace();
            int *rightCellOfFace = grid->GetRightCellOfFace();
            for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                pc[re] = pc[le];
            }
        }
        else if (bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            int *leftCellOfFace = grid->GetLeftCellOfFace();
            int *rightCellOfFace = grid->GetRightCellOfFace();
            for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                pc[re] = 0.0;
            }
        }
        else if (bcType == PHENGLEI::FARFIELD)
        {
            RDouble refflowT;
            RDouble refflowP;
            bcData->GetData("initT", &refflowT, PHDOUBLE, 1);
            bcData->GetData("initP", &refflowP, PHDOUBLE, 1);
            for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                if (FaceFlux[iFacelocal] < 0)
                {
                    pc[re] = pc[le];
                }
                else
                {
                    p[re] = refflowP;
                }
            }
        }
    }

    CommunicateAnInterfaceVar(pc);

    string GradCalcMethod = GlobalDataBase::GetStrParaFromDB("GradCalcMethod");

    GradientCalculation(grid, "pc", "dpcdx", "dpcdy", "dpcdz", GradCalcMethod);

    RDouble refPressureOfLocatePos = 0.0;
    int refPLocate = GlobalDataBase::GetIntParaFromDB("refPLocate");
    if (refPLocate == 1)
    {
        refPressureOfLocatePos = pc[0];
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        u[iCell] -= vol[iCell] / huU[iCell] * dpcdx[iCell];
        v[iCell] -= vol[iCell] / huV[iCell] * dpcdy[iCell];
        w[iCell] -= vol[iCell] / huW[iCell] * dpcdz[iCell];
        p[iCell] += PPEqRelaxCoeff * (pc[iCell] - refPressureOfLocatePos);
    }

    CommunicateAnInterfaceVar(u);
    CommunicateAnInterfaceVar(v);
    CommunicateAnInterfaceVar(w);
    CommunicateAnInterfaceVar(p);
}

void PPEqCalculator::correctBoundaryVelocity(Grid*  gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *p = reinterpret_cast<RDouble*>(grid->GetDataPtr("P"));
    RDouble *ptotal = reinterpret_cast<RDouble*>(grid->GetDataPtr("Ptotal"));
    RDouble *u = reinterpret_cast<RDouble*>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble*>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble*>(grid->GetDataPtr("W"));
    RDouble *rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
    RDouble *dpdx = (RDouble*)grid->GetDataPtr("dpdx");
    RDouble *dpdy = (RDouble*)grid->GetDataPtr("dpdy");
    RDouble *dpdz = (RDouble*)grid->GetDataPtr("dpdz");
    RDouble *dpcdx = (RDouble*)grid->GetDataPtr("dpcdx");
    RDouble *dpcdy = (RDouble*)grid->GetDataPtr("dpcdy");
    RDouble *dpcdz = (RDouble*)grid->GetDataPtr("dpcdz");
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *dunb = reinterpret_cast<RDouble*>(grid->GetDataPtr("dun"));
    RDouble *FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));
    int nTotalCell = grid->GetNTotalCell();
    RDouble *cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));
    RDouble *pc = reinterpret_cast<RDouble*>(grid->GetDataPtr("pc"));
    RDouble initRg = GlobalDataBase::GetDoubleParaFromDB("initRg");

    RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
    RDouble *huU = reinterpret_cast<RDouble *>(grid->GetDataPtr("huU"));
    RDouble *huV = reinterpret_cast<RDouble *>(grid->GetDataPtr("huV"));
    RDouble *huW = reinterpret_cast<RDouble *>(grid->GetDataPtr("huW"));
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            RDouble ub = 0.0;
            if (bcData->IsExist("initU", PHDOUBLE, 1))
            {
                bcData->GetData("initU", &ub, PHDOUBLE, 1);
            }

            RDouble vb = 0.0;
            if (bcData->IsExist("initV", PHDOUBLE, 1))
            {
                bcData->GetData("initV", &vb, PHDOUBLE, 1);
            }

            RDouble wb = 0.0;
            if (bcData->IsExist("initW", PHDOUBLE, 1))
            {
                bcData->GetData("initW", &wb, PHDOUBLE, 1);
            }

            for (std::size_t i = 0; i < faceIndex->size(); ++i)
            {
                int iFace = (*faceIndex)[i];
                int le = leftCellOfFace[iFace];
                int re = rightCellOfFace[iFace];

                u[re] = ub;
                v[re] = vb;
                w[re] = wb;
            }
        }
        else if (bcType == PHENGLEI::SYMMETRY)
        {
            for (std::size_t i = 0; i < faceIndex->size(); ++i)
            {
                int iFace = (*faceIndex)[i];
                int le = leftCellOfFace[iFace];
                int re = rightCellOfFace[iFace];

                RDouble Ud = u[le] * xfn[iFace] + v[le] * yfn[iFace] + w[le] * zfn[iFace];

                u[re] = u[le] - Ud * xfn[iFace];
                v[re] = v[le] - Ud * yfn[iFace];
                w[re] = w[le] - Ud * zfn[iFace];
            }
        }
        else if (bcType == PHENGLEI::FARFIELD)
        {
            RDouble refflowU;
            RDouble refflowV;
            RDouble refflowW;
            bcData->GetData("initU", &refflowU, PHDOUBLE, 1);
            bcData->GetData("initV", &refflowV, PHDOUBLE, 1);
            bcData->GetData("initW", &refflowW, PHDOUBLE, 1);

            for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                u[re] -= huU[le] * dpcdx[le];
                v[re] -= huV[le] * dpcdy[le];
                w[re] -= huW[le] * dpcdz[le];

                if (FaceFlux[iFacelocal] < 0)
                {
                    u[re] = refflowU;
                    v[re] = refflowV;
                    w[re] = refflowW;
                }
            }
        }
        else if (bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
        
                u[re] = u[le];
                v[re] = v[le];
                w[re] = w[le];
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET )
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                u[re] = u[le];
                v[re] = v[le];
                w[re] = w[le];
            }
        }
    }

    CommunicateAnInterfaceVar(u);
    CommunicateAnInterfaceVar(v);
    CommunicateAnInterfaceVar(w);
}

void PPEqCalculator::UpdateBCValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *p = reinterpret_cast<RDouble*>(grid->GetDataPtr("P"));
    RDouble *ptotal = reinterpret_cast<RDouble*>(grid->GetDataPtr("Ptotal"));
    RDouble *u = reinterpret_cast<RDouble*>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble*>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble*>(grid->GetDataPtr("W"));
    RDouble *rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
    RDouble *dpdx = (RDouble*)grid->GetDataPtr("dpdx");
    RDouble *dpdy = (RDouble*)grid->GetDataPtr("dpdy");
    RDouble *dpdz = (RDouble*)grid->GetDataPtr("dpdz");
    RDouble *dpcdx = (RDouble*)grid->GetDataPtr("dpcdx");
    RDouble *dpcdy = (RDouble*)grid->GetDataPtr("dpcdy");
    RDouble *dpcdz = (RDouble*)grid->GetDataPtr("dpcdz");
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *dunb = reinterpret_cast<RDouble*>(grid->GetDataPtr("dun"));
    RDouble *FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));
    int nTotalCell = grid->GetNTotalCell();
    RDouble *cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));
    RDouble *pc = reinterpret_cast<RDouble*>(grid->GetDataPtr("pc"));
    RDouble initRg = GlobalDataBase::GetDoubleParaFromDB("initRg");

    RDouble PPEqRelaxCoeff = GlobalDataBase::GetDoubleParaFromDB("PPEqRelaxCoeff");

    RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
    RDouble *InverseAp = reinterpret_cast<RDouble *>(grid->GetDataPtr("InverseAp"));
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY || bcType == PHENGLEI::INFLOW)
        {
            for (std::size_t i = 0; i < faceIndex->size(); ++i)
            {
                int iFace = (*faceIndex)[i];
                int re = rightCellOfFace[iFace];
                int le = leftCellOfFace[iFace];

                p[re] += PPEqRelaxCoeff * pc[re];
            }
        }
        else if (bcType == PHENGLEI::FARFIELD)
        {
            RDouble refMachNumber = 0.0;
            RDouble refflowU;
            RDouble refflowV;
            RDouble refflowW;
            bcData->GetData("initU", &refflowU, PHDOUBLE, 1);
            bcData->GetData("initV", &refflowV, PHDOUBLE, 1);
            bcData->GetData("initW", &refflowW, PHDOUBLE, 1);
            RDouble refflowT;
            RDouble refflowP;
            bcData->GetData("initT", &refflowT, PHDOUBLE, 1);
            bcData->GetData("initP", &refflowP, PHDOUBLE, 1);
            for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                if (FaceFlux[iFacelocal] < 0)
                {
                    p[re] = p[le];
                }
                else
                {
                    p[re] = refflowP;
                }
            }
        }
        else if (bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::PRESSURE_OUTLET)
        {
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET)
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                RDouble ptotalb = 0.0;
                if (bcData->IsExist("totalPressure", PHDOUBLE, 1))
                {
                    bcData->GetData("totalPressure", &ptotalb, PHDOUBLE, 1);
                }
                ptotal[re] = ptotalb;
                p[re] = ptotal[re] - 0.5 * rho[le] * (pow(u[re], 2) + pow(v[re], 2) + pow(w[re], 2));
            }
        }
    }

    CommunicateAnInterfaceVar(p);
}

void PPEqCalculator::solvePPEquation(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    RDouble *pc = reinterpret_cast<RDouble *> (grid->GetDataPtr("pc"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));
    RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("upperMatrixCoeff"));
    RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("lowerMatrixCoeff"));

    int solverIndex = IDX::S_IP;
    SetField(pc, 0.0, nTotalCell);

    InitializeMatrixACoeff(grid, solverIndex);

    treatBC(grid);

    constructMatrixACoeff(grid);

    constructBCoeff(grid);

    CorrectPPEb(grid);

    calculateLinearEquation(grid, IDX::S_IP, pc, diagMatrixCoeff, bCoeff, upperMatrixCoeff , lowerMatrixCoeff);

    CorrectPPEx(grid);

    correctPressureAndVelocity(grid);

    UpdateBCValue(grid);

    correctBoundaryVelocity(grid);
}

void PPEqCalculator::constructMatrixACoeff(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nIFace = grid->GetNIFace();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    RDouble *pc = reinterpret_cast<RDouble *> (grid->GetDataPtr("pc"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));
    RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("upperMatrixCoeff"));
    RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("lowerMatrixCoeff"));

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *vol = grid->GetCellVolume();
    RDouble *area = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
    RDouble *bfU = reinterpret_cast<RDouble *>(grid->GetDataPtr("bfU"));
    RDouble *bfV = reinterpret_cast<RDouble *>(grid->GetDataPtr("bfV"));
    RDouble *bfW = reinterpret_cast<RDouble *>(grid->GetDataPtr("bfW"));
    RDouble *dpdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdx"));
    RDouble *dpdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdy"));
    RDouble *dpdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdz"));
    RDouble *huU = reinterpret_cast<RDouble *>(grid->GetDataPtr("huU"));
    RDouble *huV = reinterpret_cast<RDouble *>(grid->GetDataPtr("huV"));
    RDouble *huW = reinterpret_cast<RDouble *>(grid->GetDataPtr("huW"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *mu  = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble *dun = reinterpret_cast<RDouble *>(grid->GetDataPtr("dun"));
    RDouble *rhoFace = reinterpret_cast<RDouble *>(grid->GetDataPtr("rhoFace"));
    RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));
    RDouble *ptotal = reinterpret_cast<RDouble*>(grid->GetDataPtr("Ptotal"));
    RDouble *cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        RDouble ax = area[iFace] * xfn[iFace] * (vol[le] + vol[re]) / (huU[le] + huU[re]);
        RDouble ay = area[iFace] * yfn[iFace] * (vol[le] + vol[re]) / (huV[le] + huV[re]);
        RDouble az = area[iFace] * zfn[iFace] * (vol[le] + vol[re]) / (huW[le] + huW[re]);

        RDouble dx = xcc[re] - xcc[le];
        RDouble dy = ycc[re] - ycc[le];
        RDouble dz = zcc[re] - zcc[le];

        RDouble ad = ax * dx + ay * dy + az * dz;
        RDouble a2 = ax * ax + ay * ay + az * az;

        RDouble aj = rhoFace[iFace] * a2 / ad;

        upperMatrixCoeff[iFace] = -aj;
        diagMatrixCoeff[le] += aj;

        lowerMatrixCoeff[iFace] = -aj;
        diagMatrixCoeff[re] += aj;

        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "CompressibleSIMPLE")
        {
            RDouble cRhoFace = cRho[le] * faceWeightOfLeftCell[iFace] + cRho[re] * (1.0 - faceWeightOfLeftCell[iFace]);
            RDouble speedCRho = FaceFlux[iFace] * cRhoFace / rhoFace[iFace];

            if (FaceFlux[iFace] >= 0)
            {
                diagMatrixCoeff[le] += speedCRho;
                lowerMatrixCoeff[iFace] -= speedCRho;
            }
            else
            {
                diagMatrixCoeff[re] -= speedCRho;
                upperMatrixCoeff[iFace] += speedCRho;
            }
        }
    }

    calcTransMatrixTerm(grid);
}

void PPEqCalculator::treatBC(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nIFace = grid->GetNIFace();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    RDouble *pc = reinterpret_cast<RDouble *> (grid->GetDataPtr("pc"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));
    RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("upperMatrixCoeff"));
    RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("lowerMatrixCoeff"));
    RDouble* area = grid->GetFaceArea();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *vol = grid->GetCellVolume();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    RDouble *u = reinterpret_cast<RDouble *> (grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *> (grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *> (grid->GetDataPtr("W"));
    RDouble *p = reinterpret_cast<RDouble *> (grid->GetDataPtr("P"));
    RDouble *FaceFlux = reinterpret_cast<RDouble *> (grid->GetDataPtr("FaceFlux"));
    RDouble *rho = reinterpret_cast<RDouble *> (grid->GetDataPtr("rho"));
    RDouble *huU = reinterpret_cast<RDouble *> (grid->GetDataPtr("huU"));
    RDouble *huV = reinterpret_cast<RDouble *> (grid->GetDataPtr("huV"));
    RDouble *huW = reinterpret_cast<RDouble *> (grid->GetDataPtr("huW"));
    RDouble *ptotal = reinterpret_cast<RDouble *> (grid->GetDataPtr("Ptotal"));
    RDouble *rhoFace = reinterpret_cast<RDouble*>(grid->GetDataPtr("rhoFace"));
    RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));
    RDouble *cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));
    RDouble *dpdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdx"));
    RDouble *dpdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdy"));
    RDouble *dpdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdz"));
    RDouble *bfU = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfU"));
    RDouble *bfV = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfV"));
    RDouble *bfW = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfW"));

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::FARFIELD)
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                if (FaceFlux[iFacelocal] > 0)
                {
                    RDouble ax = xfn[iFacelocal] * huU[le];
                    RDouble ay = yfn[iFacelocal] * huV[le];
                    RDouble az = zfn[iFacelocal] * huW[le];

                    RDouble dx = xfc[iFacelocal] - xcc[le];
                    RDouble dy = yfc[iFacelocal] - ycc[le];
                    RDouble dz = zfc[iFacelocal] - zcc[le];

                    RDouble ad = ax * dx + ay * dy + az * dz;
                    RDouble a2 = ax * ax + ay * ay + az * az;

                    diagMatrixCoeff[le] += rho[re] *area[iFacelocal] * a2 / ad;
                }
            }
        }
        else if (bcType == PHENGLEI::INFLOW)
        {
            string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
            if (flowSolverName == "CompressibleSIMPLE")
            {
                RDouble refflowU = 0.0;
                RDouble refflowV = 0.0;
                RDouble refflowW = 0.0;
                RDouble totalVelocity = 0.0;

                if (bcData->IsExist("initU", PHDOUBLE, 1))
                {
                    bcData->GetData("initU", &refflowU, PHDOUBLE, 1);
                }

                if (bcData->IsExist("initV", PHDOUBLE, 1))
                {
                    bcData->GetData("initV", &refflowV, PHDOUBLE, 1);
                }

                if (bcData->IsExist("initW", PHDOUBLE, 1))
                {
                    bcData->GetData("initW", &refflowW, PHDOUBLE, 1);
                }
                
                totalVelocity = sqrt(SQR(refflowU, refflowV, refflowW));

                if (bcData->IsExist("totalVelocity", PHDOUBLE, 1))
                {
                    bcData->GetData("totalVelocity", &totalVelocity, PHDOUBLE, 1);
                }

                if(totalVelocity<340.0)
                {
                    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                    {
                        int iFacelocal = (*faceIndex)[iFace];
                        int le = leftCellOfFace[iFacelocal];
                        int re = rightCellOfFace[iFacelocal];

                        if (FaceFlux[iFacelocal] >= 0)
                        {
                            diagMatrixCoeff[le] += FaceFlux[iFacelocal] / rho[le] * cRho[le];
                            lowerMatrixCoeff[iFacelocal] -= FaceFlux[iFacelocal] / rho[le] * cRho[le];
                        }
                        else
                        {
                            diagMatrixCoeff[re] -= FaceFlux[iFacelocal] / rho[re] * cRho[re];
                            upperMatrixCoeff[iFacelocal] += FaceFlux[iFacelocal] / rho[re] * cRho[re];
                        }
                    }
                }
                else
                {
                    //supersonic
                }
            }
        }
        else if (bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                RDouble ax = xfn[iFacelocal] * vol[le] / huU[le];
                RDouble ay = yfn[iFacelocal] * vol[le] / huV[le];
                RDouble az = zfn[iFacelocal] * vol[le] / huW[le];

                RDouble dx = xfc[iFacelocal] - xcc[le];
                RDouble dy = yfc[iFacelocal] - ycc[le];
                RDouble dz = zfc[iFacelocal] - zcc[le];

                RDouble ad = ax * dx + ay * dy + az * dz;
                RDouble a2 = ax * ax + ay * ay + az * az;

                diagMatrixCoeff[le] += rho[re] * area[iFacelocal] * a2 / ad;
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET)
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                RDouble ax = xfn[iFacelocal] * huU[le];
                RDouble ay = yfn[iFacelocal] * huV[le];
                RDouble az = zfn[iFacelocal] * huW[le];

                RDouble dx = xfc[iFacelocal] - xcc[le];
                RDouble dy = yfc[iFacelocal] - ycc[le];
                RDouble dz = zfc[iFacelocal] - zcc[le];

                RDouble ad = ax * dx + ay * dy + az * dz;
                RDouble a2 = ax * ax + ay * ay + az * az;

                RDouble Df= area[iFacelocal] * a2 / ad;

                if (p[re] > ptotal[re])
                {
                    p[re] = ptotal[re];
                }

                RDouble vtotal = sqrt(pow(u[re],2) + pow(v[re],2) + pow(w[re],2));
                RDouble fluxbValue = rho[re] * area[iFacelocal] * (u[re] * xfn[iFacelocal]  + v[re] * yfn[iFacelocal] + w[re] * zfn[iFacelocal]);

                diagMatrixCoeff[le] += fluxbValue * rho[le] * Df / (fluxbValue - Df * (pow(rho[le] * vtotal, 2)));
            }
        }
        else if (bcType == PHENGLEI::INTERFACE)
        {
            for (size_t iFaceLocal = 0; iFaceLocal < faceIndex->size(); ++iFaceLocal)
            {
                int iFace = (*faceIndex)[iFaceLocal];
                int le = leftCellOfFace[iFace];
                int re = rightCellOfFace[iFace];

                RDouble ax = xfn[iFace] * (vol[le] + vol[re]) / (huU[le] + huU[re]);
                RDouble ay = yfn[iFace] * (vol[le] + vol[re]) / (huV[le] + huV[re]);
                RDouble az = zfn[iFace] * (vol[le] + vol[re]) / (huW[le] + huW[re]);

                RDouble dx = xcc[re] - xcc[le];
                RDouble dy = ycc[re] - ycc[le];
                RDouble dz = zcc[re] - zcc[le];

                RDouble ad = ax * dx + ay * dy + az * dz;
                RDouble a2 = ax * ax + ay * ay + az * az;

                RDouble aj = rhoFace[iFace] * area[iFace] * a2 / ad;

                upperMatrixCoeff[iFace] = -aj;
                diagMatrixCoeff[le] += aj;

                string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
                if (flowSolverName == "CompressibleSIMPLE")
                {
                    RDouble cRhoFace = cRho[le] * faceWeightOfLeftCell[iFace] + cRho[re] * (1.0 - faceWeightOfLeftCell[iFace]);
                    RDouble speedCRho = FaceFlux[iFace] * cRhoFace / rhoFace[iFace];

                    if (FaceFlux[iFace] >= 0)
                    {
                        diagMatrixCoeff[le] += speedCRho;
                    }
                    else
                    {
                        upperMatrixCoeff[iFace] += speedCRho;
                    }
                }
            }
        }
    }
}

void PPEqCalculator::constructBCoeff(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalInFace = nTotalFace - nBoundFace;
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));

    InitializeBCoeff(grid);

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        if (bcType == PHENGLEI::FARFIELD || bcType == PHENGLEI::INFLOW || bcType == PHENGLEI::OUTFLOW
            || bcType == PHENGLEI::PRESSURE_INLET || bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                bCoeff[le] -= FaceFlux[iFacelocal];
            }
        }
        else if (bcType == PHENGLEI::INTERFACE)
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                bCoeff[le] -= FaceFlux[iFacelocal];
            }
        }
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        bCoeff[le] -= FaceFlux[iFace];
        bCoeff[re] += FaceFlux[iFace];
    }
}

void PPEqCalculator::calcTransMatrixTerm(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 0)
    {
        return;
    }

    string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
    if (flowSolverName == "CompressibleSIMPLE")
    {
        string TranCalcMethod = GlobalDataBase::GetStrParaFromDB("TranCalcMethod");
        if (TranCalcMethod == "IMPLICIT_EULER")
        {
            TransMatrixTerm_1st(grid);
        }
        else if (TranCalcMethod == "IMPLICIT_2ND_ORDER")
        {
            TransMatrixTerm_2nd(grid);
        }
    }
}

void PPEqCalculator::TransMatrixTerm_1st(Grid *gridIn) 
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *rho = reinterpret_cast<RDouble *> (grid->GetDataPtr("rho"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
    RDouble dt = GlobalDataBase::GetDoubleParaFromDB("dt");

    RDouble *cvOld = grid->GetCellVolume();
    RDouble *rhoOld = reinterpret_cast<RDouble*>(grid->GetDataPtr("rhoOld"));

    RDouble *cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));
    rhoOld = rho;
    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        diagMatrixCoeff[iCell] += cvOld[iCell] * cRho[iCell] / dt;
    }
}

void PPEqCalculator::TransMatrixTerm_2nd(Grid *gridIn) 
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *rho = reinterpret_cast<RDouble *> (grid->GetDataPtr("rho"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));
    RDouble dt = GlobalDataBase::GetDoubleParaFromDB("dt");

    RDouble *cvOld = grid->GetCellVolume();
    RDouble *rhoOld = reinterpret_cast<RDouble*>(grid->GetDataPtr("rhoOld"));
    RDouble *cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));
    RDouble *rhoOldOld = reinterpret_cast<RDouble*>(grid->GetDataPtr("rhoOldOld"));
    rhoOld = rho;
    rhoOldOld = rho;
    int timeStepNow;
    GlobalDataBase::GetData("timeStepNow", &timeStepNow, PHINT, 1);
    if (timeStepNow == 0)
    {
        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            diagMatrixCoeff[iCell] += cvOld[iCell] * cRho[iCell] / dt;
        }
    }
    else
    {
        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            diagMatrixCoeff[iCell] += 3 * cvOld[iCell] * cRho[iCell] / (2 * dt);
        }
    }
}

void PPEqCalculator::CorrectPPEx(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *pc = reinterpret_cast<RDouble *>(grid->GetDataPtr("pc"));
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    int flagOfPOutlet = 0;

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            flagOfPOutlet = 1;
            break;
        }
    }

    int totallFlagOfPOutlet = flagOfPOutlet;
    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_AllReduce(&flagOfPOutlet, &totallFlagOfPOutlet, 1, PH_SUM);
    }

    if(totallFlagOfPOutlet == 0)
    {
        RDouble xTe = 0.0;
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            xTe += pc[iCell];
        }

        RDouble globalxTe = xTe;
        RDouble GlobalTotalCells = GlobalDataBase::GetDoubleParaFromDB("GlobalTotalCells");
        if (PHMPI::GetNumberOfProcessor() > 1)
        {
            PH_AllReduce(&xTe, &globalxTe, 1, PH_SUM);
        }

        globalxTe = globalxTe / GlobalTotalCells;

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            pc[iCell] -= globalxTe;
        }
    }
}

void PPEqCalculator::CorrectPPEb(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    int flagOfPOutlet = 0;

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            flagOfPOutlet = 1;
            break;
        }
    }

    int totallFlagOfPOutlet = flagOfPOutlet;
    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_AllReduce(&flagOfPOutlet, &totallFlagOfPOutlet, 1, PH_SUM);
    }

    if(totallFlagOfPOutlet == 0)
    {
        RDouble bTe = 0.0;
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            bTe += bCoeff[iCell];
        }

        RDouble globalbTe = bTe;
        RDouble GlobalTotalCells = GlobalDataBase::GetDoubleParaFromDB("GlobalTotalCells");
        if (PHMPI::GetNumberOfProcessor() > 1)
        {
            PH_AllReduce(&bTe, &globalbTe, 1, PH_SUM);
        }

        globalbTe = globalbTe / GlobalTotalCells;

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            bCoeff[iCell] -= globalbTe;
        }
    }
}

}