#include "LineLUSGS.h"
#include "BLUSGS.h"
#include "Flux_RoeEntropyCorrection.h"
#include "Gas.h"
#include "Glb_Dimension.h"
#include "Jacobian.h"

namespace PHSPACE
{
using namespace IDX;
//! Calculate D, L, U matrixs for line implicit LUSGS methods.
void CalDLUMatrixLineLUSGS(Grid* gridIn)
{
    //! get the mesh information.
    UnstructGrid* grid = UnstructGridCast(gridIn);
   
    int nTotalCell = grid->GetNTotalCell();
    int* leftCellofFace = grid->GetLeftCellOfFace();
    int* rightCellofFace = grid->GetRightCellOfFace();

    RDouble* xfn = grid->GetFaceNormalX();
    RDouble* yfn = grid->GetFaceNormalY();
    RDouble* zfn = grid->GetFaceNormalZ();
    RDouble* area = grid->GetFaceArea();
    int* face_number_of_each_cell = grid->GetFaceNumberOfEachCell(); //Get number of faces for each cell
    int** cell2face = grid->GetCell2Face(); //Get cell to face conections

    //! get line implicit parameter, the line mesh information.
    int   nLine = grid->GetLI_nLine();
    int*  nCellsInLine = grid->GetLI_nCellsInLine();
    int** CellOfLine = grid->GetLI_CellOfLine();
    int** FaceOfLine = grid->GetLI_FaceOfLine();

    //! get the flow information.
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    bool isViscous;
    if (viscousType == INVISCID)
    {
        isViscous = false;
    }else
    {
        isViscous = true;
    }
    int nm = 5;
    RDouble** q = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));
    RDouble* dt = reinterpret_cast<RDouble*> (grid->GetDataPtr("dt"));

    //! get the three matrixs array address.
    MAT3D< RDouble >* diagMatrixMAT3D = reinterpret_cast<MAT3D< RDouble > *>(gridIn->GetDataPtr("diagMatrix"));
    MAT3D< RDouble >* lowerMatrixMAT3D = reinterpret_cast<MAT3D< RDouble > *>(gridIn->GetDataPtr("lowerMatrix"));
    MAT3D< RDouble >* upperMatrixMAT3D = reinterpret_cast<MAT3D< RDouble > *>(gridIn->GetDataPtr("upperMatrix"));

    //! create space for three matrixs. 
    if (diagMatrixMAT3D == NULL)
    {
        diagMatrixMAT3D = new MAT3D< RDouble >(nTotalCell, nm, nm);
        lowerMatrixMAT3D = new MAT3D< RDouble >(nTotalCell, nm, nm);
        upperMatrixMAT3D = new MAT3D< RDouble >(nTotalCell, nm, nm);

        grid->UpdateDataPtr("diagMatrix", diagMatrixMAT3D);
        grid->UpdateDataPtr("lowerMatrix", lowerMatrixMAT3D);
        grid->UpdateDataPtr("upperMatrix", upperMatrixMAT3D);
    }

    RDouble*** diagMatrix = diagMatrixMAT3D->GetT3D();
    RDouble*** lowerMatrix = lowerMatrixMAT3D->GetT3D();
    RDouble*** upperMatrix = upperMatrixMAT3D->GetT3D();

    //£¡ get the viscous flow parameters.
    RDouble prandtlLaminar = 0.0, prandtlTurbulence = 0.0;
    RDouble* visl = NULL, * vist = NULL;
    RDouble* normalDistanceC2C = NULL;
    if (isViscous)
    {
        visl = reinterpret_cast<RDouble*> (grid->GetDataPtr("visl"));
        vist = reinterpret_cast<RDouble*> (grid->GetDataPtr("vist"));
        prandtlLaminar = GlobalDataBase::GetDoubleParaFromDB("prl");
        prandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("prt");
        normalDistanceC2C = grid->GetNormalDistanceOfC2C();
    }

    //! define three counter for three matrixs.
    int counterDMatrix = 0;
    int counterLMatrix = 0;
    int counterUMatrix = 0;

    //! define some temporary arrays.
    RDouble* xyzfn = new RDouble[3]();
    RDouble* qThis = new RDouble[nm]();
    RDouble* qNeighbor = new RDouble[nm]();
    RDouble* qAverage = new RDouble[nm]();
    RDouble** jacobianMatrix = NewPointer2<RDouble>(nm, nm);
    RDouble** convMatrix = NewPointer2<RDouble>(nm, nm);
    RDouble** visMatrix = NewPointer2<RDouble>(nm, nm); //Jv * 2  // jacobian matrix of viscous flux

    bool D_roeMatrix, D_eigvMatrix;
    bool LU_roeMatrix, LU_eigvMatrix;
    bool D_visAllMatrix, D_visEigvMatrix;
    bool LU_visAllMatrix, LU_visEigvMatrix;

    //adjust the coefficient according the mach number.
    if (refMachNumber < 1.5)
    {
        D_roeMatrix = true;
        D_eigvMatrix = false;

        LU_roeMatrix = true;
        LU_eigvMatrix = false;

        D_visAllMatrix = true;
        D_visEigvMatrix = false;

        LU_visAllMatrix = true;
        LU_visEigvMatrix = false;
    }
    else
    {
        D_roeMatrix = false;
        D_eigvMatrix = true;

        LU_roeMatrix = false;
        LU_eigvMatrix = true;

        D_visAllMatrix = false;
        D_visEigvMatrix = true;

        LU_visAllMatrix = false;
        LU_visEigvMatrix = true;
    }

    for (int iLine = 0; iLine < nLine; ++iLine)
    {
        //! get the number of cells belong to iLine.
        int nCellsOfThisLine = nCellsInLine[iLine];

        for (int jCell = 0; jCell < nCellsOfThisLine; ++jCell)
        {
            //! get the cell number of jth cell of iLine.
            int cellNum = CellOfLine[iLine][jCell];

            RDouble diagValue = one / dt[cellNum];

            //init the diag_matrix matrix.
            for (int i = 0; i < nm; ++i)
            {
                for (int j = 0; j < nm; ++j)
                {
                    diagMatrix[counterDMatrix][i][j] = 0.0;
                }
                diagMatrix[counterDMatrix][i][i] = diagValue;
            }

            //! get primtive variables of this cell.
            qThis[IR] = q[IR][cellNum];
            qThis[IU] = q[IU][cellNum];
            qThis[IV] = q[IV][cellNum];
            qThis[IW] = q[IW][cellNum];
            qThis[IP] = q[IP][cellNum];

            for (int jFace = 0; jFace < face_number_of_each_cell[cellNum]; jFace++)
            {
                int faceNum = cell2face[cellNum][jFace];
                int cellNeighbor = leftCellofFace[faceNum] + rightCellofFace[faceNum] - cellNum;

                //! judge the face normal direction. if cellNum is owner cell, the face normal direction remain unchanged.
                //!                     
                xyzfn[0] = xfn[faceNum];
                xyzfn[1] = yfn[faceNum];
                xyzfn[2] = zfn[faceNum];
                if (leftCellofFace[faceNum] != cellNum)
                {
                    xyzfn[0] = -xfn[faceNum];
                    xyzfn[1] = -yfn[faceNum];
                    xyzfn[2] = -zfn[faceNum];
                }

                //! get primtive variables of this cell.
                qNeighbor[IR] = q[IR][cellNeighbor];
                qNeighbor[IU] = q[IU][cellNeighbor];
                qNeighbor[IV] = q[IV][cellNeighbor];
                qNeighbor[IW] = q[IW][cellNeighbor];
                qNeighbor[IP] = q[IP][cellNeighbor];

                //! calculate the face q value.
                qAverage[IR] = half * (qThis[IR] + q[IR][cellNeighbor]);
                qAverage[IU] = half * (qThis[IU] + q[IU][cellNeighbor]);
                qAverage[IV] = half * (qThis[IV] + q[IV][cellNeighbor]);
                qAverage[IW] = half * (qThis[IW] + q[IW][cellNeighbor]);
                qAverage[IP] = half * (qThis[IP] + q[IP][cellNeighbor]);

                //calculate the viscous D matrix.
                if (isViscous)
                {
                    RDouble normalDistance = normalDistanceC2C[faceNum];
                    RDouble visAverageLaminar = half * (visl[cellNum] + visl[cellNeighbor]);
                    RDouble visAverageTurbulence = half * (vist[cellNum] + vist[cellNeighbor]);

                    if (D_visAllMatrix)
                    {
                        RDouble cp = one / ((refGama - one) * refMachNumber * refMachNumber);
                        RDouble kcp = (visAverageLaminar / prandtlLaminar + visAverageTurbulence / prandtlTurbulence) * cp;
                        RDouble y0 = area[faceNum] / normalDistance / refReNumber;
                        RDouble viscosity = visAverageLaminar + visAverageTurbulence;
                        
                        CalVisJacobianMatrix(visMatrix, xyzfn, qAverage, refGama, refMachNumber, viscosity, kcp, y0);
                        for (int i = 0; i < nm; ++i) 
                        {
                            for (int j = 0; j < nm; ++j) 
                            {
                                diagMatrix[counterDMatrix][i][j] += visMatrix[i][j];
                            }
                        }
                    }

                    if (D_visEigvMatrix)
                    {
                        RDouble viscosity = visAverageLaminar + visAverageTurbulence;
                        RDouble viscousSpectrumRadius = area[faceNum] * viscosity / (qAverage[IR] * normalDistance * refReNumber + SMALL);
                        for (int m = 0; m < nm; ++m)
                        {
                            diagMatrix[counterDMatrix][m][m] += viscousSpectrumRadius;
                        }
                    }
                }
                
                //! Jacobian of convective flux zero!!!
                if (D_roeMatrix)
                {
                    CalConvJacobianMatrix(convMatrix, xyzfn, qThis, refGama, area[faceNum]);
                    CalRoeMatrix(jacobianMatrix, qThis, qNeighbor, xyzfn, refGama, area[faceNum]);
                }
                else
                {
                    for (int i = 0; i < nm; ++i)
                    {
                        for (int j = 0; j < nm; ++j)
                        {
                            jacobianMatrix[i][j] = 0.0;
                            convMatrix[i][j]     = 0.0;
                        }
                    }
                }

                if (D_eigvMatrix)
                {
                    //! Eigenvalues of convective flux      
                    RDouble rl  = qThis[IR];
                    RDouble ul  = qThis[IU];
                    RDouble vl  = qThis[IV];
                    RDouble wl  = qThis[IW];
                    RDouble pl  = qThis[IP];
                    RDouble vnl = xyzfn[0] * ul + xyzfn[1] * vl + xyzfn[2] * wl;

                    RDouble rr  = qNeighbor[IR];
                    RDouble ur  = qNeighbor[IU];
                    RDouble vr  = qNeighbor[IV];
                    RDouble wr  = qNeighbor[IW];
                    RDouble pr  = qNeighbor[IP];
                    RDouble vnr = xyzfn[0] * ur + xyzfn[1] * vr + xyzfn[2] * wr;

                    RDouble pm = half * (pl + pr);
                    RDouble rm = half * (rl + rr);
                    RDouble vn = half * (vnl + vnr);
                    RDouble cm = sqrt(refGama * pm / rm);

                    RDouble inviscidSpectrumRadius = half * area[faceNum] * (ABS(vn) + cm);
                    for (int m = 0; m < nm; ++m)
                    {
                        jacobianMatrix[m][m] += inviscidSpectrumRadius;
                    }
                }
                
                for (int i = 0; i < nm; ++i)
                {
                    for (int j = 0; j < nm; ++j)
                    {
                        diagMatrix[counterDMatrix][i][j] += convMatrix[i][j] +  jacobianMatrix[i][j];
                    }
                }
            }
            //! for upper tridiagonal matrix U.
            if (jCell < nCellsOfThisLine - 1)
            {
                int cellNeighbor = CellOfLine[iLine][jCell + 1];
                int faceNum = FaceOfLine[iLine][jCell + 1];

                //! judge the face normal direction. if cellNum is owner cell, the face normal direction remain unchanged.
                xyzfn[0] = xfn[faceNum];
                xyzfn[1] = yfn[faceNum];
                xyzfn[2] = zfn[faceNum];
                if (leftCellofFace[faceNum] != cellNum)
                {
                    xyzfn[0] = -xfn[faceNum];
                    xyzfn[1] = -yfn[faceNum];
                    xyzfn[2] = -zfn[faceNum];
                }

                qNeighbor[IR] = q[IR][cellNeighbor];
                qNeighbor[IU] = q[IU][cellNeighbor];
                qNeighbor[IV] = q[IV][cellNeighbor];
                qNeighbor[IW] = q[IW][cellNeighbor];
                qNeighbor[IP] = q[IP][cellNeighbor];

                qAverage[IR] = half * (qThis[IR] + q[IR][cellNeighbor]);
                qAverage[IU] = half * (qThis[IU] + q[IU][cellNeighbor]);
                qAverage[IV] = half * (qThis[IV] + q[IV][cellNeighbor]);
                qAverage[IW] = half * (qThis[IW] + q[IW][cellNeighbor]);
                qAverage[IP] = half * (qThis[IP] + q[IP][cellNeighbor]);

                if (LU_roeMatrix)
                {
                    CalRoeMatrix(jacobianMatrix, qThis, qNeighbor, xyzfn, refGama, area[faceNum]);
                }
                else
                {
                    for (int i = 0; i < nm; ++i)
                    {
                        for (int j = 0; j < nm; ++j)
                        {
                            jacobianMatrix[i][j] = 0.0;
                        }
                    }
                }

                if (LU_eigvMatrix)
                { 
                    //Eigenvalues of convective flux
                    RDouble rl  = qThis[IR];
                    RDouble ul  = qThis[IU];
                    RDouble vl  = qThis[IV];
                    RDouble wl  = qThis[IW];
                    RDouble pl  = qThis[IP];
                    RDouble vnl = xyzfn[0] * ul + xyzfn[1] * vl + xyzfn[2] * wl;

                    RDouble rr  = qNeighbor[IR];
                    RDouble ur  = qNeighbor[IU];
                    RDouble vr  = qNeighbor[IV];
                    RDouble wr  = qNeighbor[IW];
                    RDouble pr  = qNeighbor[IP];
                    RDouble vnr = xyzfn[0] * ur + xyzfn[1] * vr + xyzfn[2] * wr;

                    RDouble pm = half * (pl + pr);
                    RDouble rm = half * (rl + rr);
                    RDouble vn = half * (vnl + vnr);
                    RDouble cm = sqrt(refGama * pm / rm);

                    RDouble inviscidSpectrumRadius = half * area[faceNum] * (ABS(vn) + cm);

                    for (int m = 0; m < nm; ++m)
                    {
                        jacobianMatrix[m][m] += inviscidSpectrumRadius;
                    }
                }

                if (isViscous)
                {
                    RDouble dist = normalDistanceC2C[faceNum];
                    RDouble visAverageLaminar = half * (visl[cellNum] + visl[cellNeighbor]);
                    RDouble visAverageTurbulence = half * (vist[cellNum] + vist[cellNeighbor]);
                    RDouble viscosity = visAverageLaminar + visAverageTurbulence;
                    if (LU_visAllMatrix)
                    {
                        RDouble cp = one / ((refGama - one) * refMachNumber * refMachNumber);
                        RDouble kcp = (visAverageLaminar / prandtlLaminar + visAverageTurbulence / prandtlTurbulence) * cp;
                        RDouble y0 = area[faceNum] / dist / refReNumber;
                        CalVisJacobianMatrix(visMatrix, xyzfn, qAverage, refGama, refMachNumber, viscosity, kcp, y0);
                        for (int i = 0; i < nm; ++i)
                        {
                            for (int j = 0; j < nm; ++j)
                            {
                                jacobianMatrix[i][j] += visMatrix[i][j];
                            }
                        }
                    }

                    if (LU_visEigvMatrix)
                    {
                        RDouble viscousSpectrumRadius = area[faceNum] * viscosity / (qAverage[IR] * dist * refReNumber + SMALL);
                        for (int m = 0; m < nm; ++m)
                        {
                            jacobianMatrix[m][m] += viscousSpectrumRadius;
                        }
                    }
                }
                
                CalConvJacobianMatrix(convMatrix, xyzfn, qNeighbor, refGama, area[faceNum]);

                for (int i = 0; i < nm; ++i)
                {
                    for (int j = 0; j < nm; ++j)
                    {
                        upperMatrix[counterUMatrix][i][j] = convMatrix[i][j] - jacobianMatrix[i][j];
                    }
                }
                ++counterUMatrix;
            }

            if (jCell > 0)
            {  
                //for lower tridiagonal
                int cellNeighbor = CellOfLine[iLine][jCell - 1];
                int faceNum = FaceOfLine[iLine][jCell];

                //! judge the face normal direction. if cellNum is owner cell, the face normal direction remain unchanged.
                xyzfn[0] = xfn[faceNum];
                xyzfn[1] = yfn[faceNum];
                xyzfn[2] = zfn[faceNum];
                if (leftCellofFace[faceNum] != cellNum)
                {
                    xyzfn[0] = -xfn[faceNum];
                    xyzfn[1] = -yfn[faceNum];
                    xyzfn[2] = -zfn[faceNum];
                }

                qNeighbor[IR] = q[IR][cellNeighbor];
                qNeighbor[IU] = q[IU][cellNeighbor];
                qNeighbor[IV] = q[IV][cellNeighbor];
                qNeighbor[IW] = q[IW][cellNeighbor];
                qNeighbor[IP] = q[IP][cellNeighbor];

                qAverage[IR] = half * (qThis[IR] + q[IR][cellNeighbor]);
                qAverage[IU] = half * (qThis[IU] + q[IU][cellNeighbor]);
                qAverage[IV] = half * (qThis[IV] + q[IV][cellNeighbor]);
                qAverage[IW] = half * (qThis[IW] + q[IW][cellNeighbor]);
                qAverage[IP] = half * (qThis[IP] + q[IP][cellNeighbor]);

                if (LU_roeMatrix)
                {
                    
                    CalRoeMatrix(jacobianMatrix, qThis, qNeighbor, xyzfn, refGama, area[faceNum]);
                }
                else
                {
                    for (int i = 0; i < nm; ++i)
                    {
                        for (int j = 0; j < nm; ++j)
                        {
                            jacobianMatrix[i][j] = 0.0;
                        }
                    }
                }

                if (LU_eigvMatrix)
                { 

                    RDouble rl  = qThis[IR];
                    RDouble ul  = qThis[IU];
                    RDouble vl  = qThis[IV];
                    RDouble wl  = qThis[IW];
                    RDouble pl  = qThis[IP];
                    RDouble vnl = xyzfn[0] * ul + xyzfn[1] * vl + xyzfn[2] * wl;

                    RDouble rr  = qNeighbor[IR];
                    RDouble ur  = qNeighbor[IU];
                    RDouble vr  = qNeighbor[IV];
                    RDouble wr  = qNeighbor[IW];
                    RDouble pr  = qNeighbor[IP];
                    RDouble vnr = xyzfn[0] * ur + xyzfn[1] * vr + xyzfn[2] * wr;

                    RDouble pm = half * (pl + pr);
                    RDouble rm = half * (rl + rr);
                    RDouble vn = half * (vnl + vnr);
                    RDouble cm = sqrt(refGama * pm / rm);

                    RDouble inviscidSpectrumRadius = half * area[faceNum] * (ABS(vn) + cm);

                    for (int m = 0; m < nm; ++m)
                    {
                        jacobianMatrix[m][m] += inviscidSpectrumRadius;
                    }
                }

                if (isViscous)
                {
                    RDouble normalDistance = normalDistanceC2C[faceNum];
                    RDouble visAverageLaminar = half * (visl[cellNum] + visl[cellNeighbor]);
                    RDouble visAverageTurbulence = half * (vist[cellNum] + vist[cellNeighbor]);
                    RDouble viscosity = visAverageLaminar + visAverageTurbulence;
                    if (LU_visAllMatrix)
                    {
                        RDouble cp = one / ((refGama - one) * refMachNumber * refMachNumber);
                        RDouble kcp = (visAverageLaminar / prandtlLaminar + visAverageTurbulence / prandtlTurbulence) * cp;
                        RDouble y0 = area[faceNum] / normalDistance / refReNumber;
                        
                        CalVisJacobianMatrix(visMatrix, xyzfn, qAverage, refGama, refMachNumber, viscosity, kcp, y0);
                        for (int i = 0; i < nm; ++i)
                        {
                            for (int j = 0; j < nm; ++j)
                            {
                                jacobianMatrix[i][j] += visMatrix[i][j];
                            }
                        }
                    }

                    if (LU_visEigvMatrix)
                    {
                        RDouble viscousSpectrumRadius = area[faceNum] * viscosity / (qAverage[IR] * normalDistance * refReNumber + SMALL);
                        for (int m = 0; m < nm; ++m)
                        {
                            jacobianMatrix[m][m] += viscousSpectrumRadius;
                        }
                    }
                }
                
                CalConvJacobianMatrix(convMatrix, xyzfn, qNeighbor, refGama, area[faceNum]);

                for (int i = 0; i < nm; ++i)
                {
                    for (int j = 0; j < nm; ++j)
                    {
                        lowerMatrix[counterLMatrix][i][j] = convMatrix[i][j] - jacobianMatrix[i][j];
                    }
                }
                ++counterLMatrix;
            }
            ++counterDMatrix;
        }
    }

    DelPointer2(jacobianMatrix);
    DelPointer2(convMatrix);
    DelPointer2(visMatrix);
    delete []xyzfn;
    delete []qThis;
    delete []qNeighbor;
    delete []qAverage;
}

void GetDLUMatrixByLine(Grid* gridIn) 
{
    UnstructGrid* grid = UnstructGridCast(gridIn);

    MAT3D< RDouble >* diagMatrixMAT3D = reinterpret_cast<MAT3D< RDouble > *>(gridIn->GetDataPtr("diagMatrix"));
    MAT3D< RDouble >* lowerMatrixMAT3D = reinterpret_cast<MAT3D< RDouble > *>(gridIn->GetDataPtr("lowerMatrix"));
    MAT3D< RDouble >* upperMatrixMAT3D = reinterpret_cast<MAT3D< RDouble > *>(gridIn->GetDataPtr("upperMatrix"));

    RDouble*** diagMatrix = diagMatrixMAT3D->GetT3D();
    RDouble*** lowerMatrix = lowerMatrixMAT3D->GetT3D();
    RDouble*** upperMatrix = upperMatrixMAT3D->GetT3D();

    //line implicit parameter
    int   nLine = grid->GetLI_nLine();
    int* nCellsInLine = grid->GetLI_nCellsInLine();
    int nm = 5;
    RDouble**** DMatrixStartByLine = reinterpret_cast<RDouble****>(gridIn->GetDataPtr("DMatrixStartByLine"));
    RDouble**** UMatrixStartByLine = reinterpret_cast<RDouble****>(gridIn->GetDataPtr("UMatrixStartByLine"));
    RDouble**** LMatrixStartByLine = reinterpret_cast<RDouble****>(gridIn->GetDataPtr("LMatrixStartByLine"));
    if (DMatrixStartByLine == NULL)
    {
        //calculate LU of line block matrix && diag_matrix and LU
        DMatrixStartByLine = new RDouble * **[nLine];
        UMatrixStartByLine = new RDouble * **[nLine];
        LMatrixStartByLine = new RDouble * **[nLine];

        DMatrixStartByLine[0] = diagMatrix;
        UMatrixStartByLine[0] = upperMatrix;
        LMatrixStartByLine[0] = lowerMatrix - 1; //for start at c2(c1 in c++) in lower line matrix

        for (int iLine = 1; iLine < nLine; ++iLine)
        {
            DMatrixStartByLine[iLine] = DMatrixStartByLine[iLine - 1] + nCellsInLine[iLine - 1];
            UMatrixStartByLine[iLine] = UMatrixStartByLine[iLine - 1] + nCellsInLine[iLine - 1] - 1;
            LMatrixStartByLine[iLine] = LMatrixStartByLine[iLine - 1] + nCellsInLine[iLine - 1] - 1;
        }

        grid->UpdateDataPtr("DMatrixStartByLine", DMatrixStartByLine);
        grid->UpdateDataPtr("UMatrixStartByLine", UMatrixStartByLine);
        grid->UpdateDataPtr("LMatrixStartByLine", LMatrixStartByLine);
    }

    for (int iLine = 0; iLine < nLine; ++iLine)
    {
        int nCellsOfThisLine = nCellsInLine[iLine];
        BlockLUInverse<RDouble>(DMatrixStartByLine[iLine], UMatrixStartByLine[iLine], LMatrixStartByLine[iLine], nCellsOfThisLine, nm);
    }
}

void CalIsShockFace(Grid* gridIn, int* IsShockFace)
{
    UnstructGrid* grid = UnstructGridCast(gridIn);

    int nTotalFace = grid->GetNTotalFace();
    RDouble* xfn = grid->GetFaceNormalX();
    RDouble* yfn = grid->GetFaceNormalY();
    RDouble* zfn = grid->GetFaceNormalZ();
    int* leftCellofFace = grid->GetLeftCellOfFace();
    int* rightCellofFace = grid->GetRightCellOfFace();

    RDouble* primitiveVarFarfield = reinterpret_cast <RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble refPressure = primitiveVarFarfield[IP];

    RDouble** q = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));

    int le, re;
    RDouble dp, t, vx, vy, vz;
    RDouble ThdShock = 0.5;    // threshold for shock face
    for (int iFace = 0; iFace < nTotalFace; ++iFace)
    {
        IsShockFace[iFace] = 0;
    }
    for (int iFace = 0; iFace < nTotalFace; ++iFace)
    {
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        // average velocity
        vx = q[IU][le] + q[IU][re];
        vy = q[IV][le] + q[IV][re];
        vz = q[IW][le] + q[IW][re];
        t = vx * xfn[iFace] + vy * yfn[iFace] + vz * zfn[iFace];
        if (t > 0)
        {
            dp = q[IP][re] - q[IP][le];
        }
        else
        {
            dp = q[IP][le] - q[IP][re];
        }
        dp /= refPressure;
        if (dp > ThdShock)
        {
            IsShockFace[iFace] = 1;
        }
    }
}

//! Calculate viscous flux Jacobian for Block LUSGS /Line Implicit LUSGS.
void CalVisJacobianMatrix(RDouble** matrix, RDouble* xyzfn, RDouble* q, RDouble refGama, RDouble refMachNumber, RDouble vis, RDouble kcp, RDouble y0)
{
    RDouble nx = xyzfn[0];
    RDouble ny = xyzfn[1];
    RDouble nz = xyzfn[2];

    RDouble rho = q[IR];
    RDouble u   = q[IU];
    RDouble v   = q[IV];
    RDouble w   = q[IW];
    RDouble p   = q[IP];

    RDouble Vn = nx * u + ny * v + nz * w;
    RDouble halfV2 = half * (u * u + u * v + u * w);
    RDouble delta = vis * y0 / rho;

    RDouble T = refGama * refMachNumber * refMachNumber * p / rho;
    RDouble t0 = delta * kcp * T / vis;
    RDouble t1 = delta * kcp * T * rho * (refGama - 1) / p / vis;

    //! First row.
    matrix[0][0] = 0.0;
    matrix[0][1] = 0.0;
    matrix[0][2] = 0.0;
    matrix[0][3] = 0.0;
    matrix[0][4] = 0.0;

    //! Second row.
    matrix[1][0] = -delta * (u + third * Vn * nx);
    matrix[1][1] = delta * (third * nx * nx + 1.0);
    matrix[1][2] = delta * third * nx * ny;
    matrix[1][3] = delta * third * nx * nz;
    matrix[1][4] = 0;

    //! Third row.
    matrix[2][0] = -delta * (v + third * Vn * ny);
    matrix[2][1] = delta * third * ny * nx;
    matrix[2][2] = delta * (third * ny * ny + 1.0);
    matrix[2][3] = delta * third * ny * nz;
    matrix[2][4] = 0;

    //! Fourth row.
    matrix[3][0] = -delta * (w + third * Vn * nz);
    matrix[3][1] = delta * third * nz * nx;
    matrix[3][2] = delta * third * nz * ny;
    matrix[3][3] = delta * (third * nz * nz + 1.0);
    matrix[3][4] = 0;

    //! Fifth row.
    //! Count the effect of area and heat transfer here.
    //! Only when n and fn point to the same side of the face.
    matrix[4][0] = u * matrix[1][0] + v * matrix[2][0] + w * matrix[3][0] - t0 + t1 * halfV2;
    matrix[4][1] = u * matrix[1][1] + v * matrix[2][1] + w * matrix[3][1] - t1 * u;
    matrix[4][2] = u * matrix[1][2] + v * matrix[2][2] + w * matrix[3][2] - t1 * v;
    matrix[4][3] = u * matrix[1][3] + v * matrix[2][3] + w * matrix[3][3] - t1 * w;
    matrix[4][4] = t1;
}

void CalRoeMatrix(RDouble** matrix, RDouble* qL, RDouble* qR, RDouble* xyzfn,
    RDouble refGama, RDouble area)
{
    RDouble rl = qL[IR];
    RDouble ul = qL[IU];
    RDouble vl = qL[IV];
    RDouble wl = qL[IW];
    RDouble pl = qL[IP];

    RDouble rr = qR[IR];
    RDouble ur = qR[IU];
    RDouble vr = qR[IV];
    RDouble wr = qR[IW];
    RDouble pr = qR[IP];

    RDouble nx = xyzfn[0];
    RDouble ny = xyzfn[1];
    RDouble nz = xyzfn[2];
    RDouble t0 = refGama - 1.0;

    if (rl < 1.0e-5) rl = 1.0e-5;
    if (rr < 1.0e-5) rr = 1.0e-5;

    //! Note that not to miss the gama.
    RDouble hl = refGama * pl / rl / t0 + half * (ul * ul + vl * vl + wl * wl);
    RDouble hr = refGama * pr / rr / t0 + half * (ur * ur + vr * vr + wr * wr);

    RDouble ratio = sqrt(rr / rl);
    RDouble coef = 1.0 / (1.0 + ratio);

    RDouble u = (ul + ur * ratio) * coef;
    RDouble v = (vl + vr * ratio) * coef;
    RDouble w = (wl + wr * ratio) * coef;
    RDouble h = (hl + hr * ratio) * coef;
    RDouble halfV2 = half * (u * u + v * v + w * w);
    RDouble Vn = nx * u + ny * v + nz * w;

    RDouble c2 = t0 * (h - halfV2);
    RDouble c = sqrt(ABS(c2));

    RDouble eigv1 = ABS(Vn);
    RDouble eigv2 = ABS(Vn + c);
    RDouble eigv3 = ABS(Vn - c);

    //! Entropy fix.
    RDouble RoeEntropyFixCoef1, RoeEntropyFixCoef2;
    RoeEntropyFixCoef1 =0.1;
    RoeEntropyFixCoef2 =0.1;
    RDouble eigv_max = MAX(eigv2, eigv3);
    RDouble tmp0 = eigv_max * RoeEntropyFixCoef1;
    RDouble tmp1 = eigv_max * RoeEntropyFixCoef2;

    //! Count the effects of area and 0.5 in the Dissipation Fluxes here.
    eigv1 = half * area * MAX(tmp0, eigv1);
    eigv2 = half * area * MAX(tmp1, eigv2);
    eigv3 = half * area * MAX(tmp1, eigv3);

    RDouble y0 = eigv1;
    RDouble y1 = half * (2.0 * eigv1 - eigv2 - eigv3) / c2;
    RDouble y2 = half * (eigv2 - eigv3) / c;

    //! First row.
    matrix[0][0] = y0 - halfV2 * t0 * y1 - Vn * y2;
    matrix[0][1] = u * t0 * y1 + nx * y2;
    matrix[0][2] = v * t0 * y1 + ny * y2;
    matrix[0][3] = w * t0 * y1 + nz * y2;
    matrix[0][4] = -t0 * y1;

    //! Second row.
    matrix[1][0] = (nx * Vn * c2 - u * halfV2 * t0) * y1 + (nx * halfV2 * t0 - u * Vn) * y2;
    matrix[1][1] = (u * u * t0 - c2 * nx * nx) * y1 + (nx * u - nx * u * t0) * y2 + y0;
    matrix[1][2] = (v * u * t0 - c2 * ny * nx) * y1 + (ny * u - nx * v * t0) * y2;
    matrix[1][3] = (w * u * t0 - c2 * nz * nx) * y1 + (nz * u - nx * w * t0) * y2;
    matrix[1][4] = -u * t0 * y1 + nx * t0 * y2;

    //! Third row.
    matrix[2][0] = (ny * Vn * c2 - v * halfV2 * t0) * y1 + (ny * halfV2 * t0 - v * Vn) * y2;
    matrix[2][1] = (u * v * t0 - c2 * nx * ny) * y1 + (nx * v - ny * u * t0) * y2;
    matrix[2][2] = (v * v * t0 - c2 * ny * ny) * y1 + (ny * v - ny * v * t0) * y2 + y0;
    matrix[2][3] = (w * v * t0 - c2 * nz * ny) * y1 + (nz * v - ny * w * t0) * y2;
    matrix[2][4] = -v * t0 * y1 + ny * t0 * y2;

    //! Fourth row.
    matrix[3][0] = (nz * Vn * c2 - w * halfV2 * t0) * y1 + (nz * halfV2 * t0 - w * Vn) * y2;
    matrix[3][1] = (u * w * t0 - c2 * nx * nz) * y1 + (nx * w - nz * u * t0) * y2;
    matrix[3][2] = (v * w * t0 - c2 * ny * nz) * y1 + (ny * w - nz * v * t0) * y2;
    matrix[3][3] = (w * w * t0 - c2 * nz * nz) * y1 + (nz * w - nz * w * t0) * y2 + y0;
    matrix[3][4] = -w * t0 * y1 + nz * t0 * y2;

    //! Fifth row.
    matrix[4][0] = (-halfV2 * halfV2 * t0 + Vn * Vn * c2 - halfV2 * c2) * y1 + (halfV2 * Vn * (t0 - 1.0) - c2 * Vn / t0) * y2;
    matrix[4][1] = (u * c2 + u * halfV2 * t0 - nx * Vn * c2) * y1 + (nx * halfV2 + nx * c2 / t0 - u * Vn * t0) * y2;
    matrix[4][2] = (v * c2 + v * halfV2 * t0 - ny * Vn * c2) * y1 + (ny * halfV2 + ny * c2 / t0 - v * Vn * t0) * y2;
    matrix[4][3] = (w * c2 + w * halfV2 * t0 - nz * Vn * c2) * y1 + (nz * halfV2 + nz * c2 / t0 - w * Vn * t0) * y2;
    matrix[4][4] = -halfV2 * t0 * y1 + Vn * t0 * y2 + half * (eigv2 + eigv3);
}


/*******************************************************************************\
  Calculate Convective flux Jacobian
\*******************************************************************************/
void CalConvJacobianMatrix(RDouble** Matrix, RDouble* xzyfn, RDouble *q, RDouble refGama, RDouble area)
{
    RDouble rho = q[IR];
    RDouble u   = q[IU];
    RDouble v   = q[IV];
    RDouble w   = q[IW];
    RDouble p   = q[IP];
    
    RDouble nx = half * xzyfn[0] * area;
    RDouble ny = half * xzyfn[1] * area;
    RDouble nz = half * xzyfn[2] * area;

    RDouble t0 = refGama - 1.0;
    RDouble halfV2 = half * (u * u + v * v + w * w);
    RDouble Vn = nx * u + ny * v + nz * w;
    RDouble h = refGama * p / rho / t0 + halfV2;

    Matrix[0][0] = 0;
    Matrix[0][1] = nx;
    Matrix[0][2] = ny;
    Matrix[0][3] = nz;
    Matrix[0][4] = 0.0;

    Matrix[1][0] = t0 * halfV2 * nx - u * Vn;
    Matrix[1][1] = (1 - t0) * u * nx + Vn;
    Matrix[1][2] = -t0 * v * nx + u * ny;
    Matrix[1][3] = -t0 * w * nx + u * nz;
    Matrix[1][4] = t0 * nx;

    Matrix[2][0] = t0 * halfV2 * ny - v * Vn;
    Matrix[2][1] = -t0 * u * ny + v * nx;
    Matrix[2][2] = (1 - t0) * v * ny + Vn;
    Matrix[2][3] = -t0 * w * ny + v * nz;
    Matrix[2][4] = t0 * ny;

    Matrix[3][0] = t0 * halfV2 * nz - w * Vn;
    Matrix[3][1] = -t0 * u * nz + w * nx;
    Matrix[3][2] = -t0 * v * nz + w * ny;
    Matrix[3][3] = (1.0 - t0) * w * nz + Vn;
    Matrix[3][4] = t0 * nz;

    Matrix[4][0] = (t0 * halfV2 - h) * Vn;
    Matrix[4][1] = h * nx - t0 * u * Vn;
    Matrix[4][2] = h * ny - t0 * v * Vn;
    Matrix[4][3] = h * nz - t0 * w * Vn;
    Matrix[4][4] = refGama * Vn;
}

//! Calculate the Flux without eigiv value.
void FluxLUSGS3DNoEig(RDouble* flux, RDouble* qNeighbor, RDouble* DQ, RDouble* xyzfn, RDouble refGama)
{
    using namespace IDX;

    double nx, ny, nz;
    nx = xyzfn[0];
    ny = xyzfn[1];
    nz = xyzfn[2];

    double r, u, v, w, p, H, V2, Vn;
    r = qNeighbor[IR];
    u = qNeighbor[IU];
    v = qNeighbor[IV];
    w = qNeighbor[IW];
    p = qNeighbor[IP];

    V2 = u * u + v * v + w * w;
    Vn = nx * u + ny * v + nz * w;
    H = (refGama / (refGama - 1)) * (p / r) + half * V2;

    RDouble b1 = -Vn * DQ[IR] + nx * DQ[IRU] + ny * DQ[IRV] + nz * DQ[IRW];
    RDouble b2 = (refGama - 1) * (half * V2 * DQ[IR] - u * DQ[IRU] - v * DQ[IRV] - w * DQ[IRW] + DQ[IRE]);

    flux[IR] = Vn * DQ[IR] + b1;
    flux[IRU] = Vn * DQ[IRU] + u * b1 + nx * b2;
    flux[IRV] = Vn * DQ[IRV] + v * b1 + ny * b2;
    flux[IRW] = Vn * DQ[IRW] + w * b1 + nz * b2;
    flux[IRE] = Vn * DQ[IRE] + H * b1 + Vn * b2;
}

//! Calculate the |A|dQ - Roe delta Flux for Block LUSGS /Line Implicit LUSGS.
void RoeDeltaFlux(RDouble* flux, RDouble* qL, RDouble* qR, RDouble* DQ, RDouble* xyzfn, RDouble refGama, RDouble RoeEntropyFixCoef1, RDouble RoeEntropyFixCoef2)
{
    RDouble nx, ny, nz;

    RDouble r, u, v, w, h, c, Vn, halfV2, absVel, c2, ratio, coef;
    RDouble rl, ul, vl, wl, pl, hl, v2l;
    RDouble rr, ur, vr, wr, pr, hr, v2r;
    RDouble z0, z1, X1, X2, b1, b2;

    RDouble eigv1, eigv2, eigv3/*,eigv_max,tmp0,tmp1,tmp2*/;

    using namespace IDX;

    //! Get the face normal.
    nx = xyzfn[0];
    ny = xyzfn[1];
    nz = xyzfn[2];

    //! Get the the left and right cell value of the face for roe flux.
    rl = qL[IR];
    ul = qL[IU];
    vl = qL[IV];
    wl = qL[IW];
    pl = qL[IP];

    rr = qR[IR];
    ur = qR[IU];
    vr = qR[IV];
    wr = qR[IW];
    pr = qR[IP];

    if (rl <= 1.0e-5) rl = 1.0e-5;
    if (rr <= 1.0e-5) rr = 1.0e-5;

    v2l = ul * ul + vl * vl + wl * wl;
    v2r = ur * ur + vr * vr + wr * wr;

    RDouble gamm1 = refGama - 1.0;
    RDouble gaminv = refGama / gamm1;
    hl = gaminv * pl / rl + half * v2l;
    hr = gaminv * pr / rr + half * v2r;

    ratio = sqrt(rr / rl);
    coef = 1.0 / (1.0 + ratio);

    r = sqrt(rl * rr);
    u = (ul + ur * ratio) * coef;
    v = (vl + vr * ratio) * coef;
    w = (wl + wr * ratio) * coef;
    h = (hl + hr * ratio) * coef;

    halfV2 = half * (u * u + v * v + w * w);
	absVel = sqrt(two * halfV2);
    Vn = nx * u + ny * v + nz * w;

    c2 = gamm1 * (h - halfV2);
    c = sqrt(ABS(c2));

    eigv1 = ABS(Vn);
    eigv2 = ABS(Vn + c);
    eigv3 = ABS(Vn - c);

    Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, 0.0, 0.0, absVel, Vn, c, 3, RoeEntropyFixCoef1, RoeEntropyFixCoef2);

    z0 = u * DQ[IRU] + v * DQ[IRV] + w * DQ[IRW];
    z1 = nx * DQ[IRU] + ny * DQ[IRV] + nz * DQ[IRW];

    X1 = (2 * eigv1 - eigv2 - eigv3) / (2 * c2);
    X2 = (eigv2 - eigv3) / (2 * c);

    b1 = -halfV2 * DQ[IR] * gamm1 * X1 - Vn * DQ[IR] * X2 + z0 * gamm1 * X1 + z1 * X2 - DQ[IRE] * gamm1 * X1;
    b2 = Vn * c2 * DQ[IR] * X1 + halfV2 * DQ[IR] * gamm1 * X2 - c2 * z1 * X1 - z0 * gamm1 * X2 + DQ[IRE] * gamm1 * X2;

    flux[IR]  = eigv1 * DQ[IR] + b1 * 1 + b2 * 0;
    flux[IRU] = eigv1 * DQ[IRU] + b1 * u + b2 * nx;
    flux[IRV] = eigv1 * DQ[IRV] + b1 * v + b2 * ny;
    flux[IRW] = eigv1 * DQ[IRW] + b1 * w + b2 * nz;
    flux[IRE] = (eigv2 + eigv3) / 2 * DQ[IRE] + b1 * halfV2 + b2 * Vn + c2 * z0 * X1
        + c2 * z1 * X2 / gamm1 - halfV2 * c2 * DQ[IR] * X1 - c2 * Vn * DQ[IR] * X2 / gamm1;
}

//! Calculate the viscous delta Flux for Block LUSGS /Line Implicit LUSGS.
void ViscousDeltaFlux(RDouble*  flux, RDouble* q, RDouble* dQ, RDouble* xyzfn, RDouble refGama,
    RDouble vis, RDouble kcp, RDouble T)
{
    RDouble r = q[IR];
    RDouble u = q[IU];
    RDouble v = q[IV];
    RDouble w = q[IW];
    RDouble p = q[IP];

    RDouble nx = xyzfn[0];
    RDouble ny = xyzfn[1];
    RDouble nz = xyzfn[2];

    RDouble Vn = u * nx + v * ny + w * nz;
    RDouble V2 = u * u + v * v + w * w;

    RDouble t0 = nx * dQ[IRU] + ny * dQ[IRV] + nz * dQ[IRW];
    RDouble t1 = q[IU] * dQ[IRU] + q[IV] * dQ[IRV] + q[IW] * dQ[IRW];

    RDouble b1 = -dQ[IR] / r;
    RDouble b2 = 1 / r;
    RDouble b3 = (t0 - dQ[IR] * Vn) / 3 / r;

    RDouble y0 = kcp * T / vis / r;
    RDouble y1 = kcp * T * (refGama - 1) / vis / p;

    flux[IR] = 0.0;
    flux[IRU] = b1 * u + b2 * dQ[IRU] + b3 * nx;
    flux[IRV] = b1 * v + b2 * dQ[IRV] + b3 * ny;
    flux[IRW] = b1 * w + b2 * dQ[IRW] + b3 * nz;
    flux[IRE] = u * flux[IRU] + v * flux[IRV] + w * flux[IRW] + (y1 * half * V2 - y0) * dQ[IR] - y1 * t1 + y1 * dQ[IRE];
}

void SolveLineLUSGSForwardInHouse(Grid* gridIn, FieldProxy* dqProxy, FieldProxy* LUplusDQ, RDouble& sweepNormal)
{
    UnstructGrid* grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int* leftCellofFace = grid->GetLeftCellOfFace();
    int* rightCellofFace = grid->GetRightCellOfFace();
    int* face_number_of_each_cell = grid->GetFaceNumberOfEachCell();
    int** cell2face = grid->GetCell2Face();

    RDouble* xfn = grid->GetFaceNormalX();
    RDouble* yfn = grid->GetFaceNormalY();
    RDouble* zfn = grid->GetFaceNormalZ();
    RDouble* area = grid->GetFaceArea();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    bool isViscous;
    if (viscousType == INVISCID)
    {
        isViscous = false;
    }
    else
    {
        isViscous = true;
    }
    int nm = 5;

    RDouble** dq = dqProxy->GetField_UNS();
    RDouble** res = reinterpret_cast <RDouble**> (gridIn->GetDataPtr("res"));
    RDouble** deltaFlux = LUplusDQ->GetField_UNS();

    int RoeEntropyFixMethod = 3;
    RDouble RoeEntropyFixCoef1, RoeEntropyFixCoef2;
    GetEntropyFixMethodAndCoefficients(RoeEntropyFixMethod, RoeEntropyFixCoef1, RoeEntropyFixCoef2);

    RDouble**** DMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("DMatrixStartByLine"));
    RDouble**** LMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("LMatrixStartByLine"));
    RDouble**** UMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("UMatrixStartByLine"));

    // Get flow variables
    RDouble** q = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));

    RDouble cp = one / ((refGama - one) * refMachNumber * refMachNumber);
    RDouble gm2 = refGama * refMachNumber * refMachNumber;

    RDouble prandtlLaminar, prandtlTurbulence;
    RDouble* visl, * vist;
    RDouble* normalDistanceC2C;
    if (isViscous)
    {
        visl = reinterpret_cast<RDouble*> (grid->GetDataPtr("visl"));
        vist = reinterpret_cast<RDouble*> (grid->GetDataPtr("vist"));
        prandtlLaminar = GlobalDataBase::GetDoubleParaFromDB("prl");
        prandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("prt");
        normalDistanceC2C = grid->GetNormalDistanceOfC2C();
    }
    
    RDouble* fluxInviscid = new RDouble[nm];
    RDouble* fluxRoe = new RDouble[nm];
    RDouble* fluxViscous = new RDouble[nm];
    RDouble* qNeighbor = new RDouble[nm];
    RDouble* xyzfn = new RDouble[5];
    RDouble* qm = new RDouble[nm];
    RDouble* dQNeighbor = new RDouble[nm];

    bool outLineCellEigvMatrix, outLineCellVisEigvMatrix;
    if (refMachNumber < 1.5)
    {
        outLineCellEigvMatrix = false;
        outLineCellVisEigvMatrix = false;
    }
    else
    {
        outLineCellEigvMatrix = true;
        outLineCellVisEigvMatrix = true;
    }

    //line implicit parameter
    int nLine = grid->GetLI_nLine();
    int* nCellsInLines = grid->GetLI_nCellsInLine();
    int** CellOfLine = grid->GetLI_CellOfLine();
    int* LineOfCell = grid->GetLI_LineOfCell();

    // Now the Forward Sweep by line.
    for (int iLine = 0; iLine < nLine; ++iLine)
    {
        int nCellsOfThisLine = nCellsInLines[iLine];

        RDouble** dqOld = NewPointer2<RDouble>(nCellsOfThisLine, nm);
        RDouble** Ux = NewPointer2<RDouble>(nCellsOfThisLine, nm);
        RDouble** rhs0 = NewPointer2<RDouble>(nCellsOfThisLine, nm);
        RDouble** dqTry = NewPointer2<RDouble>(nCellsOfThisLine, nm);

        for (int jCell = 0; jCell < nCellsOfThisLine; ++jCell)
        {
            int cellNum = CellOfLine[iLine][jCell];
            RDouble qThis[5];

            for (int m = 0; m < nm; m++)
            {
                qThis[m] = q[m][cellNum];

                rhs0[jCell][m] = 0.0;

                dqOld[jCell][m] = dq[m][cellNum];

                Ux[jCell][m] = deltaFlux[m][cellNum];

                dq[m][cellNum] = res[m][cellNum];

                //! Then reset it to zero to store the Lower forward sweep!
                deltaFlux[m][cellNum] = 0.0;
            }

            for (int jFace = 0; jFace < face_number_of_each_cell[cellNum]; jFace++)
            {
                int faceNum = cell2face[cellNum][jFace];

                int cellNeighbor = leftCellofFace[faceNum] + rightCellofFace[faceNum] - cellNum;

                //! ghost cell
                if (cellNeighbor >= nTotalCell)
                {
                    continue;
                }
                //! c_tmp in this line or in Upper
                if (LineOfCell[cellNeighbor] > iLine)
                {
                    continue;
                }
                //! c_tmp in this line or in Upper
                if (LineOfCell[cellNeighbor] == iLine)
                {
                    continue;
                }

                // Now its neighboring cell belongs to lower triangular
                if (leftCellofFace[faceNum] == cellNum)
                {
                    xyzfn[0] = xfn[faceNum];
                    xyzfn[1] = yfn[faceNum];
                    xyzfn[2] = zfn[faceNum];
                    xyzfn[3] = area[faceNum];
                    xyzfn[4] = 0.0;
                }
                else
                {
                    xyzfn[0] = -xfn[faceNum];
                    xyzfn[1] = -yfn[faceNum];
                    xyzfn[2] = -zfn[faceNum];
                    xyzfn[3] = area[faceNum];
                    xyzfn[4] = 0.0;
                }

                for (int m = 0; m < nm; m++)
                {
                    qNeighbor[m] = q[m][cellNeighbor];
                    dQNeighbor[m] = dq[m][cellNeighbor];
                }

                // Calculate everything (I call it Flux) in lower triangular
                if (outLineCellEigvMatrix)
                {
                    //! Compute the ROE averaged face primitive variables.
                    RDouble rhoRoe, uRoe, vRoe, wRoe, pRoe;
                    if (leftCellofFace[faceNum] == cellNum)
                    {
                        RoeAveragedVariables(q[IR][cellNum], q[IU][cellNum], q[IV][cellNum], q[IW][cellNum], q[IP][cellNum], refGama,
                            q[IR][cellNeighbor], q[IU][cellNeighbor], q[IV][cellNeighbor], q[IW][cellNeighbor], q[IP][cellNeighbor], refGama, rhoRoe, uRoe, vRoe, wRoe, pRoe);
                    }
                    else
                    {
                        RoeAveragedVariables(q[IR][cellNeighbor], q[IU][cellNeighbor], q[IV][cellNeighbor], q[IW][cellNeighbor], q[IP][cellNeighbor], refGama,
                            q[IR][cellNum], q[IU][cellNum], q[IV][cellNum], q[IW][cellNum], q[IP][cellNum], refGama, rhoRoe, uRoe, vRoe, wRoe, pRoe);
                    }
                    
                    RDouble facePrimitiveVariable[5];

                    facePrimitiveVariable[IR] = rhoRoe;
                    facePrimitiveVariable[IU] = uRoe;
                    facePrimitiveVariable[IV] = vRoe;
                    facePrimitiveVariable[IW] = wRoe;
                    facePrimitiveVariable[IP] = pRoe;
                    RDouble rad = GetRadInHouse(facePrimitiveVariable, xyzfn, refGama);

                    //! The inviscid Delta Flux.
                    MXDQ_STD(qNeighbor, xyzfn[0], xyzfn[1], xyzfn[2], xyzfn[3], xyzfn[4], refGama, dQNeighbor, fluxInviscid, nm, 5, rad, -1);
                    for (int m = 0; m < nm; m++)
                    {
                        rhs0[jCell][m] += fluxInviscid[m];
                    }
                }
                else
                {
                    FluxLUSGS3DNoEig(fluxInviscid, qNeighbor, dQNeighbor, xyzfn, refGama);
                    RoeDeltaFlux(fluxRoe, qThis, qNeighbor, dQNeighbor, xyzfn, refGama, RoeEntropyFixCoef1, RoeEntropyFixCoef2);
                    for (int m = 0; m < nm; m++)
                    {
                        rhs0[jCell][m] += half * area[faceNum] * (fluxInviscid[m] - fluxRoe[m]);
                    }
                }

                if (isViscous)
                {
                    RDouble distance = normalDistanceC2C[faceNum];
                    RDouble viscousLaminarFace = half * (visl[cellNum] + visl[cellNeighbor]);
                    RDouble viscousTurbulenceFace = half * (vist[cellNum] + vist[cellNeighbor]);
                    RDouble viscosity = viscousLaminarFace + viscousTurbulenceFace;
                    if (outLineCellVisEigvMatrix)
                    {
                        RDouble density = half * (qThis[IR] + qNeighbor[IR]);
                        RDouble viscousSpectrumRadius = area[faceNum] * viscosity / (density * distance * refReNumber + SMALL);
                        for (int m = 0; m < nm; m++)
                        {
                            rhs0[jCell][m] -= viscousSpectrumRadius * dQNeighbor[m];
                        }
                    }
                    else
                    {
                        RDouble kcp = (viscousLaminarFace / prandtlLaminar + viscousTurbulenceFace / prandtlTurbulence) * cp;
                        
                        for (int m = 0; m < nm; m++)
                        {
                            qm[m] = half * (qThis[m] + qNeighbor[m]);
                        }
                        RDouble visSpecRad = viscosity * area[faceNum] / (distance * refReNumber);
                        RDouble tm = gm2 * qm[IP] / qm[IR];
                        ViscousDeltaFlux(fluxViscous, qm, dQNeighbor, xyzfn, refGama, viscosity, kcp, tm);

                        for (int m = 0; m < nm; m++)
                        {
                            rhs0[jCell][m] -= visSpecRad * fluxViscous[m];
                        }
                    }
                }
            }

            for (int m = 0; m < nm; m++)
            {
                dq[m][cellNum] = dq[m][cellNum] - Ux[jCell][m] - rhs0[jCell][m];
                dqTry[jCell][m] = dq[m][cellNum];
                deltaFlux[m][cellNum] += rhs0[jCell][m];
            }
        }

        BlockLUSolve<RDouble>(dqTry, DMatrixStartByLine[iLine], UMatrixStartByLine[iLine], LMatrixStartByLine[iLine], nCellsOfThisLine, nm);

        for (int jCell = 0; jCell < nCellsOfThisLine; ++jCell)
        {
            int cellNum = CellOfLine[iLine][jCell];
            for (int m = 0; m < nm; m++)
            {
                dq[m][cellNum] = dqTry[jCell][m];
                sweepNormal += SQR(dq[m][cellNum] - dqOld[jCell][m]);
            }
        }

        DelPointer2(dqOld);
        DelPointer2(dqTry);
        DelPointer2(Ux);
        DelPointer2(rhs0);
    }

    SetGhostDQLUSGSInHouse(grid, dq, q);

    delete [] qNeighbor;
    delete [] dQNeighbor;
    delete [] xyzfn;
    delete [] qm;
    delete [] fluxInviscid;
    delete [] fluxViscous;
    delete [] fluxRoe;
}

void SolveLineLUSGSBackwardInHouse(Grid* gridIn, FieldProxy* dqProxy, FieldProxy* LUplusDQ, RDouble& sweepNormal)
{
    UnstructGrid* grid = UnstructGridCast(gridIn);
    int* leftCellofFace = grid->GetLeftCellOfFace();
    int* rightCellofFace = grid->GetRightCellOfFace();
    int* face_number_of_each_cell = grid->GetFaceNumberOfEachCell();
    int** cell2face = grid->GetCell2Face();

    RDouble* xfn = grid->GetFaceNormalX();
    RDouble* yfn = grid->GetFaceNormalY();
    RDouble* zfn = grid->GetFaceNormalZ();
    RDouble* area = grid->GetFaceArea();

    RDouble** dq = dqProxy->GetField_UNS();
    RDouble** res = reinterpret_cast <RDouble**> (grid->GetDataPtr("res"));
    RDouble** deltaFlux = LUplusDQ->GetField_UNS();

    //! get the flow information.
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    bool isViscous;
    if (viscousType == INVISCID)
    {
        isViscous = false;
    }else
    {
        isViscous = true;
    }
    int nm = 5;
    int RoeEntropyFixMethod = 3;
    RDouble RoeEntropyFixCoef1, RoeEntropyFixCoef2;
    GetEntropyFixMethodAndCoefficients(RoeEntropyFixMethod, RoeEntropyFixCoef1, RoeEntropyFixCoef2);

    bool outLineCellEigvMatrix, outLineCellVisEigvMatrix;
    if (refMachNumber < 1.5)
    {
        outLineCellEigvMatrix = false;
        outLineCellVisEigvMatrix = false;
    }
    else
    {
        outLineCellEigvMatrix = true;
        outLineCellVisEigvMatrix = true;
    }

    RDouble**** DMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("DMatrixStartByLine"));
    RDouble**** LMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("LMatrixStartByLine"));
    RDouble**** UMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("UMatrixStartByLine"));

    //! Get flow variables
    RDouble** q = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));

    RDouble cp = one / ((refGama - one) * refMachNumber * refMachNumber);
    RDouble gm2 = refGama * refMachNumber * refMachNumber;

    RDouble prandtlLaminar, prandtlTurbulence;
    RDouble* visl, * vist;
    RDouble* normalDistanceC2C;
    if (isViscous)
    {
        visl = reinterpret_cast<RDouble*> (grid->GetDataPtr("visl"));
        vist = reinterpret_cast<RDouble*> (grid->GetDataPtr("vist"));
        prandtlLaminar = GlobalDataBase::GetDoubleParaFromDB("prl");
        prandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("prt");
        normalDistanceC2C = grid->GetNormalDistanceOfC2C();
    }

    //line implicit parameter
    int nLine = grid->GetLI_nLine();
    int* nCellsInLines = grid->GetLI_nCellsInLine();
    int** CellOfLine = grid->GetLI_CellOfLine();
    int* LineOfCell = grid->GetLI_LineOfCell();
    
    RDouble* fluxInviscid = new RDouble[nm];
    RDouble* fluxRoe = new RDouble[nm];
    RDouble* fluxViscous = new RDouble[nm];
    RDouble* qNeighbor = new RDouble[nm];
    RDouble* xyzfn = new RDouble[5];
    RDouble* qm = new RDouble[nm];
    RDouble* dQNeighbor = new RDouble[nm];
    RDouble* qThis = new RDouble[nm];

    //! Backward Sweep
    for (int iLine = nLine - 1; iLine >= 0; --iLine)
    {
        int nCellsOfThisLine = nCellsInLines[iLine];
        RDouble** dqOld = NewPointer2<RDouble>(nCellsOfThisLine, nm);
        RDouble** Lx = NewPointer2<RDouble>(nCellsOfThisLine, nm);
        RDouble** rhs0 = NewPointer2<RDouble>(nCellsOfThisLine, nm);
        RDouble** dqTry = NewPointer2<RDouble>(nCellsOfThisLine, nm);

        for (int jCell = nCellsOfThisLine - 1; jCell >= 0; --jCell)
        {
            int cellNum = CellOfLine[iLine][jCell];

            for (int m = 0; m < nm; m++)
            {
                qThis[m] = q[m][cellNum];

                dqOld[jCell][m] = dq[m][cellNum];

                Lx[jCell][m] = deltaFlux[m][cellNum];

                rhs0[jCell][m] = 0.0;

                dq[m][cellNum] = res[m][cellNum];
                
                deltaFlux[m][cellNum] = 0.0;
            }

            for (int j = 0; j < face_number_of_each_cell[cellNum]; j++)
            {
                int faceNum = cell2face[cellNum][j];

                int cellNeighbor = leftCellofFace[faceNum] + rightCellofFace[faceNum] - cellNum;

                //! c_tmp in this line or in lower
                if (LineOfCell[cellNeighbor] < iLine)
                {
                    continue;
                }
                //! c_tmp in this line or in lower
                if (LineOfCell[cellNeighbor] == iLine)
                {
                    continue;
                }

                //! Now its neighboring cell belongs to lower triangular
                xyzfn[0] = xfn[faceNum];
                xyzfn[1] = yfn[faceNum];
                xyzfn[2] = zfn[faceNum];
                xyzfn[3] = area[faceNum];
                xyzfn[4] = 0.0;
                if (leftCellofFace[faceNum] != cellNum)
                {
                    xyzfn[0] = -xfn[faceNum];
                    xyzfn[1] = -yfn[faceNum];
                    xyzfn[2] = -zfn[faceNum];
                    xyzfn[3] = area[faceNum];
                    xyzfn[4] = 0.0;
                }

                for (int m = 0; m < nm; m++)
                {
                    qNeighbor[m] = q[m][cellNeighbor];
                    dQNeighbor[m] = dq[m][cellNeighbor];
                }

                //! Calculate everything (I call it Flux) in lower triangular
                if (outLineCellEigvMatrix)
                {
                    //! Compute the ROE averaged face primitive variables.
                    RDouble rhoRoe, uRoe, vRoe, wRoe, pRoe;
                    RoeAveragedVariables(q[IR][cellNum], q[IU][cellNum], q[IV][cellNum], q[IW][cellNum], q[IP][cellNum], refGama,
                        q[IR][cellNeighbor], q[IU][cellNeighbor], q[IV][cellNeighbor], q[IW][cellNeighbor], q[IP][cellNeighbor], refGama, rhoRoe, uRoe, vRoe, wRoe, pRoe);

                    RDouble facePrimitiveVariable[5];

                    facePrimitiveVariable[IR] = rhoRoe;
                    facePrimitiveVariable[IU] = uRoe;
                    facePrimitiveVariable[IV] = vRoe;
                    facePrimitiveVariable[IW] = wRoe;
                    facePrimitiveVariable[IP] = pRoe;
                    RDouble rad = GetRadInHouse(facePrimitiveVariable, xyzfn, refGama);

                    //! The inviscid Delta Flux.
                    MXDQ_STD(qNeighbor, xyzfn[0], xyzfn[1], xyzfn[2], xyzfn[3], xyzfn[4], refGama, dQNeighbor, fluxInviscid, nm, 5, rad, -1);
                    for (int m = 0; m < nm; m++)
                    {
                        rhs0[jCell][m] += fluxInviscid[m];
                    }
                }
                else
                {
                    FluxLUSGS3DNoEig(fluxInviscid, qNeighbor, dQNeighbor, xyzfn, refGama);
                    RoeDeltaFlux(fluxRoe, qThis, qNeighbor, dQNeighbor, xyzfn, refGama, RoeEntropyFixCoef1, RoeEntropyFixCoef2);

                    for (int m = 0; m < nm; m++)
                    {
                        rhs0[jCell][m] += half * area[faceNum] * (fluxInviscid[m] - fluxRoe[m]);
                    }
                }

                if (isViscous)
                {
                    RDouble distance = normalDistanceC2C[faceNum];
                    RDouble viscousLaminarFace = half * (visl[cellNum] + visl[cellNeighbor]);
                    RDouble viscousTurbulenceFace = half * (vist[cellNum] + vist[cellNeighbor]);
                    RDouble viscosity = viscousLaminarFace + viscousTurbulenceFace;
                    if (outLineCellVisEigvMatrix)
                    {
                        RDouble density = half * (qThis[IR] + qNeighbor[IR]);
                        RDouble viscousSpectrumRadius = area[faceNum] * viscosity / (density * distance * refReNumber + SMALL);
                        for (int m = 0; m < nm; m++)
                        {
                            rhs0[jCell][m] -= viscousSpectrumRadius * dQNeighbor[m];
                        }
                    }
                    else
                    {
                        RDouble kcp = (viscousLaminarFace / prandtlLaminar + viscousTurbulenceFace / prandtlTurbulence) * cp;
                        
                        for (int m = 0; m < nm; m++)
                        {
                            qm[m] = half * (qThis[m] + qNeighbor[m]);
                        }

                        RDouble visSpecRad = viscosity * area[faceNum] / (distance * refReNumber);
                        RDouble tm = gm2 * qm[IP] / qm[IR];

                        ViscousDeltaFlux(fluxViscous, qm, dQNeighbor, xyzfn, refGama, viscosity, kcp, tm);

                        for (int m = 0; m < nm; m++)
                        {
                            rhs0[jCell][m] -= visSpecRad * fluxViscous[m];
                        }
                    }
                }
            }

            for (int m = 0; m < nm; m++)
            {
                dq[m][cellNum] = dq[m][cellNum] - Lx[jCell][m] - rhs0[jCell][m];
                deltaFlux[m][cellNum] += rhs0[jCell][m];
                dqTry[jCell][m] = dq[m][cellNum];
            }
        }

        BlockLUSolve<RDouble>(dqTry, DMatrixStartByLine[iLine], UMatrixStartByLine[iLine], LMatrixStartByLine[iLine], nCellsOfThisLine, nm);

        for (int iCell = 0; iCell < nCellsOfThisLine; ++iCell)
        {
            int cellNum = CellOfLine[iLine][iCell];
            for (int m = 0; m < nm; m++)
            {
                dq[m][cellNum] = dqTry[iCell][m];
                sweepNormal += SQR(dq[m][cellNum] - dqOld[iCell][m]);
            }
        }

        DelPointer2(dqOld);
        DelPointer2(dqTry);
        DelPointer2(rhs0);
        DelPointer2(Lx);
    }

    delete [] fluxInviscid;
    delete [] fluxRoe;
    delete [] fluxViscous;
    delete [] qNeighbor;
    delete [] xyzfn;
    delete [] qm;
    delete [] dQNeighbor;
    delete [] qThis;
}

RDouble GetRadInHouse(RDouble *prim, RDouble *nxyz, RDouble gama)
{
    using namespace IDX;

    RDouble &xfn  = nxyz[0];
    RDouble &yfn  = nxyz[1];
    RDouble &zfn  = nxyz[2];
    RDouble &area = nxyz[3];
    RDouble &vgn  = nxyz[4];

    RDouble &rm = prim[IR];
    RDouble &um = prim[IU];
    RDouble &vm = prim[IV];
    RDouble &wm = prim[IW];
    RDouble &pm = prim[IP];

    RDouble vn = xfn * um + yfn * vm + zfn * wm - vgn;
    RDouble eigen = ABS(vn) + sqrt(gama * pm / rm);
    return eigen * area;
}

}