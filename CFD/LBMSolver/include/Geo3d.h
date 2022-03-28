#pragma once

#include "LBMGeometry.h"

#include <math.h>

inline int FindCellIndex(
  const double &loc, 
  const double &origin, 
  const double &dx
){ return int(floor((loc-origin)/dx)); }

namespace Geo2d {
  inline double dotproductVV(const double *p1, const double *p2) {
    return p1[0] * p2[0] + p1[1] * p2[1];
  }
}

namespace Geo3d {
  std::pair<double, nnwlbm::geo> CalIntersection_segment_coordsSquare(
    double *intersection, // return the intersection
    const double *p1, // first point of segment
    const double *p2, // second point of segment
    const size_t di, //  face direction
    const double *lbOfFace, // left bottom of the face
    const double *rtOfFace, // right top of the face
    const double tolerance // right top of the face
  );

  double get_intersection_coords(
    const double *loc, 
    const double *norm, 
    const int &di, 
    const double *coordLine
  );

  inline void crossproduct(double *ans, const double *p1, const double *p2) {
    ans[0] = p1[1] * p2[2] - p1[2] * p2[1];
    ans[1] = p1[2] * p2[0] - p1[0] * p2[2];
    ans[2] = p1[0] * p2[1] - p1[1] * p2[0];
  }
  inline double dotproductVV(const double *p1, const double *p2) {
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
  }
  inline double triangle_area(const double *x1, const double *x2, const double *x3) {	// checked
    double v21[3], v31[3], vec[3];
    v21[0] = x2[0] - x1[0];
    v21[1] = x2[1] - x1[1];
    v21[2] = x2[2] - x1[2];
    v31[0] = x3[0] - x1[0];
    v31[1] = x3[1] - x1[1];
    v31[2] = x3[2] - x1[2];
    Geo3d::crossproduct(vec, v21, v31);
    return 0.5 * sqrt(Geo3d::dotproductVV(vec, vec));
  }
  double normalize(double *vector, const double &err = 1e-10);
  // one triangle splits into four subtriangle.
  void SplitTriangle(double *subTri, const double *triIn);
  bool triangle_3d_to_2d(double *tri2d, const double *tri3d);
  bool triangle_3d_to_2d(double *tri2d, double *matrix, const double *tri3d);
  bool on_triangle(const double *pos, const double *tri3d, const double &tolerance);
  int TriangleIntersection(double *ans, const double *triangle, 
    const double *p1, const double *p2, double tolerance);
  inline bool is_positive_side(const double *triPos, const double *triNorm, const double *pos, const double &tolerance = 1e-10) {
    double tempV[3];
    tempV[0] = pos[0] - triPos[0];
    tempV[1] = pos[1] - triPos[1];
    tempV[2] = pos[2] - triPos[2];
    return Geo3d::dotproductVV(tempV, triNorm) > tolerance; 
  }
}

