#pragma once

#include "CommonAPI.h"
#include "LBMGeometry.h"
#include "Range.h"
#include <unordered_map>
#include <set>

namespace nnwlbm {
class COMMONAPI StlGeo : public Geometry
{
private:
  std::string fileName;
  std::vector<double> triangleList;
  MaskRange<int, 3> colorMap;
  Box boundingBox;
  Range<3> nodeRange;
  RangeIndex<3> gridIndex;
  double cellDx;
  std::unordered_multimap<size_t, size_t> cellToTri;
  std::unordered_multimap<size_t, size_t> triToCell;
  int numOfConnectedDomain;

  inline double *get_triangle(const size_t &index) { return triangleList.data() + 12 * index; }
  bool GenConnectedDomain();
public:
  virtual std::string name() const { return "StlGeo"; }
  explicit StlGeo() : Geometry(3), boundingBox(3), numOfConnectedDomain(0),
    gridIndex(Range<3>(), RangeIndex<3>::inputNodeRange) {}
  virtual int CalRelativeLocation(
    const double*, 
    const double &tolerance = 1e-9
  ) const;
  virtual int CalRelativeLocationOrigin(
    const double *, 
    const double &tolerance = 1e-9
  ) const;
  virtual bool cal_normal(
    double *norm, 
    const double *pos, 
    const double &tolerance = 1e-9
  ) const;
  virtual bool CalIntersection(
    double *intersection, 
    const double *p1, 
    const double *p2, 
    const double &tolerance = 1e-9
  ) const;
  bool ReadFile(
    const std::string &fname, 
    const double &scaleFactor = 1.0
  );
  inline void set_triangleVec(std::vector<double> &triVec) { triangleList = triVec; }
  inline size_t get_triangleNum() const { return triangleList.size() / 12; }
  inline const double *GetTrianglePos(const size_t &index) const { return triangleList.data() + 12 * index + 3; }
  inline const double *GetTriangleNorm(const size_t &index) const { return triangleList.data() + 12 * index; }
  void GetTriIdList(
    std::set<size_t> &triIdList, 
    const size_t &cellId
  ) const;
  bool GenBoundingBox(const double &minL);
  const Box &GetBoundingBox() const { return boundingBox; }
  void WriteColorMap() const;
  std::shared_ptr<MaskRange<int,3>> sharpen_boundary() const;
};

}
