//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      MeshDeformationRBF.cpp (Win32 ONLY!!!)
//! @brief     This file gives an example of  the function of points deformation by RBF method.
//!            The serial points deformation procedure is illustrated, and if you want the parallel version, 
//!            you may contact the following author to get supports. 
//!            For this serial version of RBF points deformation function, user only need simple input, that is:\n
//!                in put: 1. All of the original points (x, y, z) that have not been deformed;\n
//!                        2. The deformed points (x, y, z) and them deformation vector (dx, dy, dz), 
//!                           which are usually the deformed boundary points such as solid wall deformation.\n
//!                        3. Several controlling parameters, such as the influence radius.
//!                out put: All of the newly deformed points.
//! @author    Bell.

//#include "PHengLEI.h"
#include "ParallelRBF.h"
#include <fstream>

using namespace PHSPACE;

void ReadPoints(PointVector *pointsCoord);
void ReadDeformedPoints(PointVector *deformedPoints, PointVector *deformation);
void DumpNewPoints(PointVector *pointsCoord);

int main()
{
    //InitPHengLEI();

    //*********************************************************************************************
    //             Start preparing the input of the RBF module.
    //*********************************************************************************************
    //! dimension: 2 or 3.
    //! The RBF influence radius parameter: multiply it to the maximum boundary deformation distance
    //! is the influence radius used in RBF method (influenceRadius in following).
    //! the points that out of the radius will not been deformed.
    //! 10.0-50.0 is recommended.
    int dimension = 2;
    double influenceRadius = 100.0;

    ParallelRBF * parallelRBF = new ParallelRBF(dimension, influenceRadius); //! Create the parallel RBF method.
    
    //! Number of controlling points, that are used for basis of the deformation.
    //! 5-10 points is recommended.
    int numberOfReferenceCP = 5;
    parallelRBF->SetNumberOfReferenceCP(numberOfReferenceCP);

    //! Number of zones in this processors and total global zones. For serial computing, all set to be '1'.
    int nLocalZones = 1, nGlobalZones = 1;
    parallelRBF->SetNumberOfZone(nLocalZones, nGlobalZones);

    //! Read all of the original points (x, y, z) that have not been deformed from file.
    //! The first dimension is the number of zones (for serial it is '1');
    //! The second dimension is the number of points index in the corresponding zone.
    //! These points, usually the whole field, will be deformed by the following 'deformedPoints' && 'deformation'.
    PointVector pointsCoordinates;
    pointsCoordinates.resize(1);
    ReadPoints(&pointsCoordinates);
    parallelRBF->SetPointsCoordinates(&pointsCoordinates);

    //! Read the deformed points (x, y, z) and them deformation vector (dx, dy, dz), 
    //! which are usually the deformed boundary points such as solid wall deformation.\n
    //! The first dimension is the number of zones (for serial it is '1');
    //! The second dimension is the number of points index in the corresponding zone.
    //! These points, will used to deform the 'pointsCoordinates'.
    PointVector deformedPoints, deformation;
    deformedPoints.resize(1);
    deformation.resize(1);
    ReadDeformedPoints(&deformedPoints, &deformation);
    parallelRBF->SetBoundaryDeformedPoints(&deformedPoints);
    parallelRBF->SetDeformedVector(&deformation);

    //! Here, we assume each of the 'pointsCoordinates' need to deform.
    int nTotalPoints = pointsCoordinates[0].size();
    vector< vector<bool> > needDeform;
    needDeform.resize(1);
    needDeform[0].resize(nTotalPoints);
    for (int iNode = 0; iNode < nTotalPoints; ++ iNode)
    {
        needDeform[0][iNode] = true;
    }
    parallelRBF->SetIsPointsNeedDeform(&needDeform);  
    
    //*********************************************************************************************
    //                          The kernel part: RBF run
    //*********************************************************************************************
    //! RBF mesh deformation begin.     
    parallelRBF->Run();
    //! RBF mesh deformation End


    //*********************************************************************************************
    //             Check the deformed points after deformation
    //*********************************************************************************************
    //! Dump out the new deformed points.
    DumpNewPoints(&pointsCoordinates);

   

    //! Exit PHengLEI.
    //ExitPHengLEI();
    return 0; 
}

void ReadPoints(PointVector *pointsCoord)
{
    int nTotalNode;
    ifstream infile("../in/pointsCoordinates.dat", ios::in);
    infile >> nTotalNode;
    (*pointsCoord)[0].resize(nTotalNode);
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        double x0, y0, z0;
        infile >> x0 >> y0 >> z0;

        (*pointsCoord)[0][iNode].SetX(x0);
        (*pointsCoord)[0][iNode].SetY(y0);
        (*pointsCoord)[0][iNode].SetZ(z0);
    }
    infile.close();
}

void ReadDeformedPoints(PointVector *deformedPoints, PointVector *deformation)
{
    int nDeformedNode;
    ifstream infile("../in/pointsDeformation.dat", ios::in);
    infile >> nDeformedNode;
    (*deformedPoints)[0].resize(nDeformedNode);
    (*deformation)[0].resize(nDeformedNode);    
    
    for (int iNode = 0; iNode < nDeformedNode; ++ iNode)
    {
        double x0, y0, z0, dx, dy, dz;
        infile >> x0 >> y0 >> z0 >> dx >> dy >> dz;

        (*deformedPoints)[0][iNode].SetX(x0);
        (*deformedPoints)[0][iNode].SetY(y0);
        (*deformedPoints)[0][iNode].SetZ(z0);

        (*deformation)[0][iNode].SetX(dx);
        (*deformation)[0][iNode].SetY(dy);
        (*deformation)[0][iNode].SetZ(dz);
    }
    infile.close();
}

void DumpNewPoints(PointVector *pointsCoord)
{
    int nTotalNode = (*pointsCoord)[0].size();
    ofstream outfile("../out/pointsCoordinatesNew.dat", ios::out);
    outfile.setf(ios::scientific);
    outfile.precision(8);
    //outfile << nTotalNode << endl;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        double x0 = (*pointsCoord)[0][iNode].X();
        double y0 = (*pointsCoord)[0][iNode].Y();
        double z0 = (*pointsCoord)[0][iNode].Z();

        outfile << x0 << " " << y0 << " " << z0 << endl;
    }
    outfile.close();
}