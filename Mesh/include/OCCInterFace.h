#pragma once
#ifdef USE_OCC
#include <iostream>
#include <IGESControl_Reader.hxx>
#include <TColStd_HSequenceOfTransient.hxx>
#include <TopoDS.hxx>
#include <Standard_Macro.hxx>
#include <Handle_IGESDefs_MacroDef.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepTools.hxx>
#include <TopExp.hxx>
#include <Standard.hxx>
#include <Standard_Transient_proto.hxx>
#include <IGESToBRep_TopoCurve.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include "Precision.h"
using namespace std;
using namespace PHSPACE;

class OCCInterFace
{
private:
    // Dimension of the problem: 2 or 3.
    int dimension;

    // the shape
    TopoDS_Shape shape;

    // all the (sub) TopoDS_Shapes in 'shape'
    TopTools_IndexedMapOfShape faceMap, edgeMap, vertexMap, solidMap, shellMap, wireMap;

private:
    void Scale(RDouble scaling);
    void buildLists();
    void addShapeToLists(TopoDS_Shape & shape);
    void projectionTest();

public:
    OCCInterFace();
    ~OCCInterFace();

    void loadIGES(const char * fileName, int geometryUnit);
    void ComputeNearestPointOnGeometry(const gp_Pnt & pointIn, gp_Pnt & pointNearest, 
        RDouble & projectDistance, bool & isProjected, RDouble * borderTolerancePercent);
    void SetDimension(int dimIn) { this->dimension = dimIn; }
    TopTools_IndexedMapOfShape & GetFaceMap() { return this->faceMap; }
    void PrintIGESBox();
};

#endif