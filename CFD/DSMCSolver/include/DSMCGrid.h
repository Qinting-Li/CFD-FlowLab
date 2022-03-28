//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      DSMCGrid.h
//! @brief     Grid definition for DSMC.
//! @author    Li Jin.

#pragma once
#include "Geo_UnstructGrid.h"

namespace PHSPACE
{

class DSMCUnstructGrid
{
private:
    int nMolecule,moleculeStartID;
    int nWallFace,nInletFace,nOutletFace,nSymmetryFace,nInnerFace;
    int *boundaryFaceMap,*wallFaceMap,*inletFaceMap,*outletFaceMap,*symmetryFaceMap,*innerFaceMap;
    int *boundaryFaceBlock;

    //int *nOutVect_WL,*nOutVect_FR,*nOutVect_SM;
    int *subCellLocation;
    int *wallFace2BoundaryFace, *inletFace2BoundaryFace, *outletFace2BoundaryFace, *symmetryFace2BoundaryFace, *interFace2BoundaryFace;
    RDouble * cellWeight;
public:
    DSMCUnstructGrid();
    ~DSMCUnstructGrid();

    void SetnMolecule(int nMolecule){this->nMolecule = nMolecule;}
    int GetnMolecule() const {return this->nMolecule ;}

    void CalCellWeight(UnstructGrid * grid,RDouble);
    void SetCellWeight(RDouble * cellWeight){this->cellWeight = cellWeight;}
    RDouble * GetCellWeight() const {return this->cellWeight ;}

    void SetnWallFace(int size){this->nWallFace = size;}
    int GetnWallFace() const {return this->nWallFace ;}

    void SetnInletFace(int nInletFace){this->nInletFace = nInletFace;}
    int GetnInletFace() const {return this->nInletFace ;}

    void SetnOutletFace(int nOutletFace){this->nOutletFace = nOutletFace;}
    int GetnOutletFace() const {return this->nOutletFace ;}

    void SetnSymmetryFace(int nSymmetryFace){this->nSymmetryFace = nSymmetryFace;}
    int GetnSymmetryFace() const {return this->nSymmetryFace ;}

    void SetnInnerFace(int nInnerFace){this->nInnerFace = nInnerFace;}
    int GetnInnerFace() const {return this->nInnerFace ;}

    void BoundaryCount(UnstructGrid * grid);

    void SetWallFace2BoundaryFace(int *wallFace2BoundaryFace){this->wallFace2BoundaryFace = wallFace2BoundaryFace; }
    int * GetWallFace2BoundaryFace(){return this->wallFace2BoundaryFace;}
    
    void SetInletFace2BoundaryFace(int *inletFace2BoundaryFace){this->inletFace2BoundaryFace = inletFace2BoundaryFace; }
    int * GetInletFace2BoundaryFace(){return this->inletFace2BoundaryFace;}

    void SetOutletFace2BoundaryFace(int *outletFace2BoundaryFace){this->outletFace2BoundaryFace = outletFace2BoundaryFace; }
    int * GetOutletFace2BoundaryFace(){return this->outletFace2BoundaryFace;}

    void SetSymmetryFace2BoundaryFace(int *symmetryFaceMap){this->symmetryFaceMap = symmetryFaceMap; }
    int * GetSymmetryFace2BoundaryFace(){return this->symmetryFaceMap;}

    void SetInterFace2BoundaryFace(int *interFace2BoundaryFace){this->interFace2BoundaryFace = interFace2BoundaryFace; }
    int * GetInterFace2BoundaryFace(){return this->interFace2BoundaryFace;}
};

    RDouble CalMaxCellVolume();
    RDouble CalMinCellVolume();
    void CalReferenceCellVolume();
    RDouble GetReferenceCellVolume();
}



