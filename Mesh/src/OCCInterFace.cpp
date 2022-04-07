#ifdef USE_OCC
#include <OCCInterFace.h>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Vertex.hxx>
#include <Geom_Curve.hxx>
#include <BRep_Tool.hxx>
#include <gp_Pnt.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <ShapeAnalysis.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include "Math_BasisFunction.h"
#include "TK_Log.h"
#include "PHHeader.h"
#include "TK_Time.h"
OCCInterFace::OCCInterFace()
{
    dimension = 3;
}

OCCInterFace::~OCCInterFace()
{

}

void OCCInterFace::loadIGES(const char * fileName, int geometryUnit)
{
    PHSPACE::PrintToWindow("Start loading IGES ...", "\n");

    IGESControl_Reader reader;
    reader.ReadFile((char*) fileName);
    reader.NbRootsForTransfer();
    reader.TransferRoots();
    shape = reader.OneShape();
    BRepTools::Clean(shape);

//  reader.PrintTransferInfo(IFSelect_FailOnly, IFSelect_ResultCount);
//  reader.PrintTransferInfo(IFSelect_FailOnly, IFSelect_Mapping);

    //! Geometry unit.\n
    //!   1 -- meter;\n
    //!   2 -- millimeter;\n
    //!   3 -- inch.
    if(geometryUnit == 1)
    {
        Scale(1.0/25.4);
    }else if(geometryUnit == 2)
    {
        Scale(1.0/25400);
    }else if(geometryUnit == 3)
    {

    }
    

    buildLists();

    PHSPACE::PrintToWindow("End loading IGES ...", "\n");
}

void OCCInterFace::Scale(RDouble scaling)
{
    gp_Trsf transfer;
    transfer.SetScaleFactor(scaling);
    BRepBuilderAPI_Transform transform(this->shape, transfer);
    this->shape = transform.Shape();
}

void OCCInterFace::buildLists()
{
    solidMap.Clear();
    shellMap.Clear();
    faceMap.Clear();
    wireMap.Clear();
    edgeMap.Clear();
    vertexMap.Clear();

    addShapeToLists(shape);
}

void OCCInterFace::addShapeToLists(TopoDS_Shape & shape)
{
    // Solids
    TopExp_Explorer expSolid, expShell, expFace, expWire, expEdge, expVertex;
    for(expSolid.Init(shape, TopAbs_SOLID); expSolid.More(); expSolid.Next())
    {
        TopoDS_Solid solid = TopoDS::Solid(expSolid.Current());
        if(solidMap.FindIndex(solid) < 1)
        {
            solidMap.Add(solid);

            for(expShell.Init(solid, TopAbs_SHELL); expShell.More(); expShell.Next())
            {
                TopoDS_Shell shell = TopoDS::Shell(expShell.Current());
                if(shellMap.FindIndex(shell) < 1)
                {
                    shellMap.Add(shell);

                    for(expFace.Init(shell, TopAbs_FACE); expFace.More(); expFace.Next())
                    {
                        TopoDS_Face face = TopoDS::Face(expFace.Current());
                        if(faceMap.FindIndex(face) < 1)
                        {
                            faceMap.Add(face);

                            for(expWire.Init(expFace.Current(), TopAbs_WIRE); expWire.More(); expWire.Next())
                            {
                                TopoDS_Wire wire = TopoDS::Wire(expWire.Current());
                                if(wireMap.FindIndex(wire) < 1)
                                {
                                    wireMap.Add(wire);

                                    for(expEdge.Init(expWire.Current(), TopAbs_EDGE); expEdge.More(); expEdge.Next())
                                    {
                                        TopoDS_Edge edge = TopoDS::Edge(expEdge.Current());
                                        if(edgeMap.FindIndex(edge) < 1)
                                        {
                                            edgeMap.Add(edge);

                                            for(expVertex.Init(expWire.Current(), TopAbs_VERTEX); expVertex.More(); expVertex.Next())
                                            {
                                                TopoDS_Vertex vertex = TopoDS::Vertex(expVertex.Current());
                                                if(vertexMap.FindIndex(vertex) < 1)
                                                { vertexMap.Add(vertex); }
                                            } // vertex
                                        }
                                    }// edge
                                }
                            } // wire
                        }
                    }// face
                }
            }// Shell
        } // solid
    }

    // TopExp_Explorer expSolid, expShell, expFace, expWire, expEdge, expVertex;
    // Free Shells
    for(expShell.Init(expSolid.Current(), TopAbs_SHELL, TopAbs_SOLID); expShell.More(); expShell.Next())
    {
        TopoDS_Shape shell = expShell.Current();
        if(shellMap.FindIndex(shell) < 1)
        {
            shellMap.Add(shell);

            for(expFace.Init(shell, TopAbs_FACE); expFace.More(); expFace.Next())
            {
                TopoDS_Face face = TopoDS::Face(expFace.Current());
                if(faceMap.FindIndex(face) < 1)
                {
                    faceMap.Add(face);

                    for(expWire.Init(expFace.Current(), TopAbs_WIRE); expWire.More(); expWire.Next())
                    {
                        TopoDS_Wire wire = TopoDS::Wire(expWire.Current());
                        if(wireMap.FindIndex(wire) < 1)
                        {
                            wireMap.Add(wire);

                            for(expEdge.Init(expFace.Current(), TopAbs_EDGE); expEdge.More(); expEdge.Next())
                            {
                                TopoDS_Edge edge = TopoDS::Edge(expEdge.Current());
                                if(edgeMap.FindIndex(edge) < 1)
                                {
                                    edgeMap.Add(edge);

                                    for(expVertex.Init(expEdge.Current(), TopAbs_VERTEX); expVertex.More(); expVertex.Next())
                                    {
                                        TopoDS_Vertex vertex = TopoDS::Vertex(expVertex.Current());
                                        if(vertexMap.FindIndex(vertex) < 1)
                                        { vertexMap.Add(vertex); }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // TopExp_Explorer expSolid, expShell, expFace, expWire, expEdge, expVertex;
    // Free Faces
    for(expFace.Init(shape, TopAbs_FACE, TopAbs_SHELL); expFace.More(); expFace.Next())
    {
        TopoDS_Face face = TopoDS::Face(expFace.Current());
        if(faceMap.FindIndex(face) < 1)
        {
            faceMap.Add(face);

            for(expWire.Init(expFace.Current(), TopAbs_WIRE); expWire.More(); expWire.Next())
            {
                TopoDS_Wire wire = TopoDS::Wire(expWire.Current());
                if(wireMap.FindIndex(wire) < 1)
                {
                    wireMap.Add(wire);

                    for(expEdge.Init(expWire.Current(), TopAbs_EDGE); expEdge.More(); expEdge.Next())
                    {
                        TopoDS_Edge edge = TopoDS::Edge(expEdge.Current());
                        if(edgeMap.FindIndex(edge) < 1)
                        {
                            edgeMap.Add(edge);

                            for(expVertex.Init(expEdge.Current(), TopAbs_VERTEX); expVertex.More(); expVertex.Next())
                            {
                                TopoDS_Vertex vertex = TopoDS::Vertex(expVertex.Current());
                                if(vertexMap.FindIndex(vertex) < 1)
                                { vertexMap.Add(vertex); }
                            }
                        }
                    }
                }
            }
        }
    }

    // TopExp_Explorer expSolid, expShell, expFace, expWire, expEdge, expVertex;
    // Free Wires
    for(expWire.Init(shape, TopAbs_WIRE, TopAbs_FACE); expWire.More(); expWire.Next())
    {
        TopoDS_Wire wire = TopoDS::Wire(expWire.Current());
        if(wireMap.FindIndex(wire) < 1)
        {
            wireMap.Add(wire);

            for(expEdge.Init(expWire.Current(), TopAbs_EDGE); expEdge.More(); expEdge.Next())
            {
                TopoDS_Edge edge = TopoDS::Edge(expEdge.Current());
                if(edgeMap.FindIndex(edge) < 1)
                {
                    edgeMap.Add(edge);

                    for(expVertex.Init(expEdge.Current(), TopAbs_VERTEX); expVertex.More(); expVertex.Next())
                    {
                        TopoDS_Vertex vertex = TopoDS::Vertex(expVertex.Current());
                        if(vertexMap.FindIndex(vertex) < 1)
                        { vertexMap.Add(vertex); }
                    }
                }
            }
        }
    }

    // TopExp_Explorer expSolid, expShell, expFace, expWire, expEdge, expVertex;
    // Free Edges
    for(expEdge.Init(shape, TopAbs_EDGE, TopAbs_WIRE); expEdge.More(); expEdge.Next())
    {
        TopoDS_Edge edge = TopoDS::Edge(expEdge.Current());
        if(edgeMap.FindIndex(edge) < 1)
        {
            edgeMap.Add(edge);

            for(expVertex.Init(expEdge.Current(), TopAbs_VERTEX); expVertex.More(); expVertex.Next())
            {
                TopoDS_Vertex vertex = TopoDS::Vertex(expVertex.Current());
                if(vertexMap.FindIndex(vertex) < 1)
                { vertexMap.Add(vertex); }
            }
        }
    }

    // TopExp_Explorer expSolid, expShell, expFace, expWire, expEdge, expVertex;
    // Free Vertices
    for(expVertex.Init(shape, TopAbs_VERTEX, TopAbs_EDGE); expVertex.More(); expVertex.Next())
    {
        TopoDS_Vertex vertex = TopoDS::Vertex(expVertex.Current());
        if(vertexMap.FindIndex(vertex) < 1)
        { vertexMap.Add(vertex); }
    }

}

void OCCInterFace::projectionTest()
{
    for(int iVertex = 1; iVertex <= vertexMap.Extent();  ++ iVertex)
    {
        TopoDS_Vertex vertex = TopoDS::Vertex(vertexMap(iVertex));

        gp_Pnt v0 = BRep_Tool::Pnt(vertex);
    
        cout << v0.X() << " " << v0.Y() << endl;
    }

    for (int iFace = 1; iFace <= faceMap.Extent(); ++ iFace)
    {
        TopoDS_Face face = TopoDS::Face(faceMap(iFace));
        Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
        if (!surface)
        {
            continue;
        }
        else
        {
            gp_Pnt point = surface->Value(0.5, 0.5);
            cout << "point on surface: " << point.X() <<  " " << point.Y() << " " << point.Z() << endl;
        }
    }

    for(int iEdge = 1; iEdge <= edgeMap.Extent();  ++ iEdge)
    {
        TopoDS_Edge edge = TopoDS::Edge(edgeMap(iEdge));

        Standard_Real first, last;
        Handle(Geom_Curve) curve;
        curve = BRep_Tool::Curve(edge, first, last);

        if(!curve)
        {
            continue;
        }

        gp_Pnt point(0.0, 0.5, 0.0);
        GeomAPI_ProjectPointOnCurve projectionCurve(point, curve, first, last);

        gp_Pnt projectedPoint = projectionCurve.NearestPoint();
        cout << "Nearest point: " << projectedPoint.X() << " " << projectedPoint.Y() << " " << projectedPoint.Z() << endl;
    }
}

void OCCInterFace::PrintIGESBox()
{
    if(this->dimension == 3)
    {
        RDouble xmin = PHSPACE::LARGE, xmax = - PHSPACE::LARGE;
        RDouble ymin = PHSPACE::LARGE, ymax = - PHSPACE::LARGE;
        RDouble zmin = PHSPACE::LARGE, zmax = - PHSPACE::LARGE;
        int numerOfFaceMap = faceMap.Extent();
        for (int iFace = 1; iFace <= numerOfFaceMap; ++ iFace)
        {
            TopoDS_Shape shape = faceMap(iFace);
            Bnd_Box box;
            BRepBndLib::Add(shape, box);

            RDouble x0, x1, y0, y1, z0, z1;
            box.Get(x0, y0, z0, x1, y1, z1);

            xmin = PHSPACE::MIN(xmin, x0);
            ymin = PHSPACE::MIN(ymin, y0);
            zmin = PHSPACE::MIN(zmin, z0);

            xmax = PHSPACE::MAX(xmax, x1);
            ymax = PHSPACE::MAX(ymax, y1);
            zmax = PHSPACE::MAX(zmax, z1);
        }

        ostringstream oss;
        oss << "Number of surface shape in IGES file: " << numerOfFaceMap << endl;
        oss << "Boundary box of IGES file: " << endl;
        oss << "	xmin && xmax: " << xmin << ", " << xmax << endl;
        oss << "	ymin && ymax: " << ymin << ", " << ymax << endl;
        oss << "	zmin && zmax: " << zmin << ", " << zmax << endl;
        PHSPACE::WriteLogFile(oss);
        PHSPACE::PrintToWindow(oss);
    }
}

void OCCInterFace::ComputeNearestPointOnGeometry(const gp_Pnt & pointIn, gp_Pnt & pointNearest, 
                                                 RDouble & projectDistance, bool & isProjected, RDouble * borderTolerancePercent)
{
    RDouble dist;
    isProjected = false;

    // Curves.
    if(this->dimension == 2)
    {
        int numberOfEdgeMap = edgeMap.Extent();
        for(int iEdge = 1; iEdge <= numberOfEdgeMap;  ++ iEdge)
        {
            TopoDS_Edge edge = TopoDS::Edge(edgeMap(iEdge));

            Standard_Real first, last;
            Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, first, last);
            if(!curve)
            {
                continue;
            }

            GeomAPI_ProjectPointOnCurve projectionCurve(pointIn, curve, first, last);
            if(!projectionCurve.NbPoints())
            {
                continue;
            }

            RDouble parameter = projectionCurve.LowerDistanceParameter();
            if(parameter < first || parameter > last)
            {
                continue;
            }

            gp_Pnt projectedPoint = projectionCurve.NearestPoint();
            isProjected = true;

            dist = projectionCurve.LowerDistance();

            if(dist < projectDistance)
            {
                projectDistance = dist;
                pointNearest.SetX(projectedPoint.X());
                pointNearest.SetY(projectedPoint.Y());
                pointNearest.SetZ(projectedPoint.Z());
            }
        }
    }


    // Faces.
    if(this->dimension == 3)
    {
        int numerOfFaceMap = faceMap.Extent();
        int count = 0;
        for (int iFace = 1; iFace <= numerOfFaceMap; ++ iFace)
        {
            // Step1: check if the point is in the box.
            TopoDS_Shape shape = faceMap(iFace);
            Bnd_Box box;
            BRepBndLib::Add(shape, box);

            //! Important: it is not necessary to enlarge the box, because when the point on
            //! the edge or on a plan, it need not to project it!.
            /*
            RDouble enlargeTolerence = box.SquareExtent() * 0.001;
            box.Enlarge(enlargeTolerence);
            */

            if(box.IsOut(pointIn))
            {
                ++ count;
                continue;
            }
            
            // Step2: Projection.
            TopoDS_Face face = TopoDS::Face(faceMap(iFace));
            Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
            if (!surface)
            {
                continue;
            }

            RDouble uMin = 0.0, uMax = 1.0, vMin = 0.0, vMax = 1.0;
            ShapeAnalysis::GetFaceUVBounds(face, uMin, uMax, vMin, vMax);
            GeomAPI_ProjectPointOnSurf projectionSurface(pointIn, surface, uMin, uMax, vMin, vMax);
            if(!projectionSurface.NbPoints())
            {
                continue;
            }

            RDouble initialGuess[2] = {0, 0};
            RDouble pp[2] = {initialGuess[0],initialGuess[1]};
            projectionSurface.LowerDistanceParameters(pp[0], pp[1]);

            RDouble uDelta = uMax - uMin;
            RDouble vDelta = vMax - vMin;
            RDouble uBorderToleranceMin = borderTolerancePercent[0] * uDelta;
            RDouble uBorderToleranceMax = borderTolerancePercent[1] * uDelta;
            RDouble vBorderToleranceMin = borderTolerancePercent[2] * vDelta;
            RDouble vBorderToleranceMax = borderTolerancePercent[3] * vDelta;

            if((pp[0] < uMin + uBorderToleranceMin || pp[0] > uMax - uBorderToleranceMax) || 
               (pp[1] < vMin + vBorderToleranceMin || pp[1] > vMax - vBorderToleranceMax))
            {
                continue;
            }

            //if(fabs(pp[0] - uMin) <= tolerance || fabs(pp[0] - uMax) <= tolerance) continue;
            //if(fabs(pp[1] - uMin) <= tolerance || fabs(pp[1] - uMax) <= tolerance) continue;

            gp_Pnt projectedPoint = projectionSurface.NearestPoint();
            isProjected = true;

            dist = projectionSurface.LowerDistance();

            if(dist < projectDistance)
            {
                projectDistance = dist;
                pointNearest.SetX(projectedPoint.X());
                pointNearest.SetY(projectedPoint.Y());
                pointNearest.SetZ(projectedPoint.Z());
            }
        }

        /*
        if(!isProjected)
        {
            cout << "Error: could not project point: (" 
                 << pointIn.X() << ", " << pointIn.Y() << ", " << pointIn.Z() << "), out of "
                 << count << " face in Total face: " << numerOfFaceMap
                 << endl;
            exit(0);
        }
        */
    }
}

#endif