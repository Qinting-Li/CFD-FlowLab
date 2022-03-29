
#include "CommonLibDll.h"
#include "StlGeo.h"
#include "Geo3d.h"
#include "UniformGrid.h"
#include "LBMGeometry.h"
#include <fstream>
#include <queue>
#include <array>

namespace Geo3d
{
  double get_intersection_coords(
    const double *anchor, 
    const double *faceNorm, 
    const int &di, 
    const double *coordLine
  ) {
    int di1 = (di + 1) % 3;
    int di2 = (di + 2) % 3;
    if(fabs(faceNorm[di]) > 1e-13) {
      return -(faceNorm[di1] * (coordLine[di1] - anchor[di1]) + faceNorm[di2] * (coordLine[di2] - anchor[di2]))
        /faceNorm[di] + anchor[di];
    }
    else return 1e13;
  }

  std::pair<double, nnwlbm::geo> CalIntersection_segment_coordsSquare(
    double *intersection,
    const double *p1, // first point of segment
    const double *p2, // second point of segment
    const size_t di, //  face direction
    const double *lbOfFace, // left bottom of the face
    const double *rtOfFace, // right top of the face
    const double tolerance // right top of the face
  ) {
    if(lbOfFace[di] != rtOfFace[di]) return std::make_pair(1e99, nnwlbm::geo::ndef);
    const size_t di1 = (di + 1) % 3;
    const size_t di2 = (di + 2) % 3;
    intersection[di] = lbOfFace[di];
    double unitVec[3]; 
    unitVec[0] = p2[di] - p1[di];
    unitVec[1] = p2[di1] - p1[di1];
    unitVec[2] = p2[di2] - p1[di2];
    //double length = sqrt(unitVec[0]*unitVec[0]+unitVec[1]*unitVec[1]+unitVec[2]*unitVec[2]);
    double length = normalize(unitVec, 1e-10);
    if(fabs(unitVec[0]) / length < 1e-15) {
      return std::make_pair(1e99, nnwlbm::geo::outside);
    }
    double ratio = (lbOfFace[di] - p1[di]) / unitVec[0];
    intersection[di1] = unitVec[1] * ratio + p1[di1];
    intersection[di2] = unitVec[2] * ratio + p1[di2];
    if(intersection[di1] > lbOfFace[di1] + tolerance && intersection[di1] < rtOfFace[di1] - tolerance
      && intersection[di2] > lbOfFace[di2] + tolerance && intersection[di2] < rtOfFace[di2] - tolerance) {
      return std::make_pair(ratio, nnwlbm::geo::inside);
    }
    if(intersection[di1] >= lbOfFace[di1] - tolerance && intersection[di1] <= rtOfFace[di1] + tolerance
      && intersection[di2] >= lbOfFace[di2] - tolerance && intersection[di2] <= rtOfFace[di2] + tolerance) {
      return std::make_pair(ratio, nnwlbm::geo::boundary);
    }
    return std::make_pair(ratio, nnwlbm::geo::outside);
  }

  bool is_intersected(
    const double *p1, // the first point of a line segment
    const double *p2, // the second piont of a line segment
    const double *box, // box
    const double &tolerance
  ) {
    std::pair<double, nnwlbm::geo> status[6];
    double intersection[3];
    double lbAndRt[6];
    int counterInside = 0;
    int counterBoundary = 0;
    int counterOutside = 0;
    int counterOnLine = 0;
    int counterOutsideLinePositive = 0;
    int counterOutsideLineNegtive = 0;
    for(size_t iDim = 0; iDim < 3; iDim++) {
      memcpy(lbAndRt, box, 6*sizeof(double));
      lbAndRt[iDim] = lbAndRt[iDim + 3] = box[iDim];
      status[0 + 2 * iDim] = CalIntersection_segment_coordsSquare(intersection, p1, p2, iDim, lbAndRt, lbAndRt + 3, tolerance);
      lbAndRt[iDim] = lbAndRt[iDim + 3] = box[iDim + 3];
      status[1 + 2 * iDim] = CalIntersection_segment_coordsSquare(intersection, p1, p2, iDim, lbAndRt, lbAndRt + 3, tolerance);
    }
    for(size_t i = 0; i < 6; i++) {
      if(status[i].second == nnwlbm::geo::outside)  counterOutside++;
      else if(status[i].second == nnwlbm::geo::inside) {
        counterInside++;
        if     (status[i].first < 0.0) counterOutsideLineNegtive++;
        else if(status[i].first > 1.0) counterOutsideLinePositive++;
        else                           counterOnLine++;
      }
      else if(status[i].second == nnwlbm::geo::boundary) {
        counterBoundary++;
        if     (status[i].first < 0.0) counterOutsideLineNegtive++;
        else if(status[i].first > 1.0) counterOutsideLinePositive++;
        else                           counterOnLine++;
      }
    }
    if(counterOnLine > 0) return true;
    if(counterOutsideLineNegtive > 0 && counterOutsideLinePositive > 0) return true;
    return false;
  }

  //! the maximum length of triangle must smaller than the minimum length of box
  bool is_intersected(
    const double *tri, // 
    const double *box, // box
    const double &tolerance
  ) {
    for(int iLine = 0; iLine < 3; iLine++) {
      if(is_intersected(tri + 3 * iLine, tri + 3 * ((iLine + 1) % 3), box, tolerance)) return true;
    }
    return false;
  }
  
  void SplitTriangle(double *subTri, const double *triIn) {
    double *middle[3];
    middle[0] = subTri + 27;  middle[1] = subTri + 30;  middle[2] = subTri + 33;
    for(int iDim = 0; iDim < 3; iDim++) {
      middle[0][iDim] = 0.5 * (triIn[iDim + 3 * 0] + triIn[iDim + 3 * 1]);
      middle[1][iDim] = 0.5 * (triIn[iDim + 3 * 1] + triIn[iDim + 3 * 2]);
      middle[2][iDim] = 0.5 * (triIn[iDim + 3 * 2] + triIn[iDim + 3 * 0]);
    }
    for(int iTri = 0; iTri < 3; iTri++) {
      int iTri1 = (iTri - 1 + 3) % 3;
      for(int iDim = 0; iDim < 3; iDim++) {
        subTri[iDim + 3 * (0 + 3 * iTri)] = triIn[iDim + 3 * iTri];
        subTri[iDim + 3 * (1 + 3 * iTri)] = middle[iTri][iDim];
        subTri[iDim + 3 * (2 + 3 * iTri)] = middle[iTri1][iDim];
      }
    }
  }

  double normalize(double *vec, const double &err) {
    double length = sqrt(dotproductVV(vec,vec));
    if(length < err) {
      vec[0] = vec[1] = vec[2] = 0.0;
      return 0.0;
    }
    else {
      double temp = 1.0 / length;
      vec[0] *= temp;
      vec[1] *= temp;
      vec[2] *= temp;
      return length;
    }
  }

  bool triangle_3d_to_2d(double *tri2d, double *matrix, const double *tri3d) {
    double v1[3],v2[3];
    for(int iDim = 0; iDim < 3; iDim++) {
      v1[iDim] = tri3d[iDim + 3 * 1] - tri3d[iDim + 3 * 0];
      v2[iDim] = tri3d[iDim + 3 * 2] - tri3d[iDim + 3 * 0];
    }
    crossproduct(matrix, v1, v2);
    double length = normalize(matrix);
    if(length < 1e-20) return false;
    tri2d[0 + 2 * 0] = tri2d[ 1 + 2 * 0] = 0.0;
    tri2d[0 + 2 * 1] = normalize(v1);
    tri2d[1 + 2 * 1] = 0.0;
    memcpy(matrix + 3, v1, 3 * sizeof(double));
    crossproduct(matrix + 6, matrix, v1);
    tri2d[0 + 2 * 2] = dotproductVV(v2,matrix + 3);
    tri2d[1 + 2 * 2] = dotproductVV(v2,matrix + 6);
    return true;
  }

  bool triangle_3d_to_2d(double *tri2d, const double *tri3d) {
    double matrix[3 * 3];
    return triangle_3d_to_2d(tri2d, matrix, tri3d);
  }

  bool on_triangle(
    const double *pos, 
    const double *tri3d, 
    const double &tolerance
  ) {
    double area[4];
    area[0] = triangle_area(tri3d + 0, tri3d + 3, tri3d + 6);
    if(area[0] < tolerance * tolerance) return false;
    area[1] = triangle_area(pos, tri3d + 0, tri3d + 3);
    area[2] = triangle_area(pos, tri3d + 3, tri3d + 6);
    area[3] = triangle_area(pos, tri3d + 6, tri3d + 0);
    if(fabs(area[1] + area[2] + area[3] - area[0]) > tolerance * area[0]) return false;
    else return true;
  }

  int TriangleIntersection(
    double *ans, 
    const double *triangle, 
    const double *p1, 
    const double *p2, 
    double tolerance
  ) {
    //! calculate the intersection first
    double matrix[3 * 3], tri2d[2 * 3], v1[3], v2[3], v12d[2], v22d[2];
    double area; // area1, area2, area3;
    double intersection[2];
    triangle_3d_to_2d(tri2d, matrix, triangle);
    double maxL = fabs(tri2d[0 + 2 * 1]) + sqrt(Geo2d::dotproductVV(tri2d + 4,tri2d + 4)); 
    if(maxL < tolerance) return nnwlbm::geo::outside;
    //! tolerance *= maxL;
    area = fabs(0.5 * tri2d[0 + 2 * 1] * tri2d[1 + 2 * 2]);
    v1[0] = p1[0] - triangle[0];
    v1[1] = p1[1] - triangle[1];
    v1[2] = p1[2] - triangle[2];
    v12d[0] = dotproductVV(v1, matrix + 3);
    v12d[1] = dotproductVV(v1, matrix + 6);
    double l1 = dotproductVV(v1, matrix);
    v2[0] = p2[0] - triangle[0];
    v2[1] = p2[1] - triangle[1];
    v2[2] = p2[2] - triangle[2];
    v22d[0] = dotproductVV(v2, matrix + 3);
    v22d[1] = dotproductVV(v2, matrix + 6);
    double l2 = dotproductVV(v2, matrix);
    if(fabs(l1) < tolerance) {
      if(fabs(l2) > tolerance) {
        //! p1 is on the surface
        intersection[0] = v12d[0];
        intersection[1] = v12d[1];
      }
      else {
        //! parallel or both on the surface. Whatever it is, the intersection is not well define, therefore return "outside"
        return nnwlbm::geo::outside;
      }
    }
    else {		// fabs(l1) > tolerance
      if(fabs(l2) < tolerance) {
        //! p2 is on the surface
        intersection[0] = v22d[0];
        intersection[1] = v22d[1];
      }
      else {  // fabs(l1) > tolerance, fabs(l2) > tolerance
        if(l1 * l2 > tolerance) return nnwlbm::geo::outside;	// p1 and p2 locate at the same side of the triangle.
        else {	// l1*l2 < 0.0, p1 and p2 reside at different sides of the triangle.
          l1 = fabs(l1);
          l2 = fabs(l2);
          double factor = l1 / (l1 + l2);
          intersection[0] = (1.0 - factor) * v12d[0] + (factor) * v22d[0];
          intersection[1] = (1.0 - factor) * v12d[1] + (factor) * v22d[1];
        }
      }
    }

    //                    x2
    //                    .
    //                   /|` .                     total area = a1+a2+a3
    //                  / \   ` .
    //                 /   |     ` .
    //                /    \    a2  ` .
    //               /  a3  |          ` .
    //              /    , '   ~   -      ` .  
    //         x0  /__'_____a1_________~___-_`_. x1

    for(int iDim = 0; iDim < 3; iDim++) {
      ans[iDim] = intersection[0] * matrix[iDim + 3 * 1] + intersection[1] * matrix[iDim + 3 * 2] + triangle[iDim];
    }
    if(!on_triangle(ans, triangle, tolerance)) return nnwlbm::geo::outside;
    //area1 = fabs(0.5*tri2d[0+2*1]*intersection[1]);
    //area2 = Geo2d::triangle_area(intersection, tri2d+2, tri2d+4);
    //area3 = 0.5*fabs(Geo2d::crossproduct(intersection, tri2d+4));
    ////tolerance *= maxL;
    //if(fabs(area-area1-area2-area3) > tolerance*area) return nnwlbm::geo::outside;
    else {
      if(fabs(l1) < tolerance || fabs(l2) < tolerance) return nnwlbm::geo::boundary;
      else return nnwlbm::geo::inside;
    }
  }
}

//! StlGeo
namespace nnwlbm {
bool StlGeo::ReadFile(
  const std::string &fname, 
  const double &scaleFactor
) {
  std::ifstream stlFile(fname.c_str(), std::ios::in | std::ios::binary);
  char headInfo[81];
  //! read 80 byte header
  if (stlFile) {
    fileName = fname;
    stlFile.read(headInfo, 80);
    headInfo[80] = '\0';
    std::cout << "STL file comment " << headInfo << std::endl;
    char nTriRaw[4];
    stlFile.read(nTriRaw, 4);
    unsigned numTri = *((unsigned*)nTriRaw);
    triangleList.resize(12 * numTri);
    for(unsigned i = 0; i < numTri; i++) {
      char triRaw[50];
      stlFile.read(triRaw, 50);
      double *ptr = triangleList.data() + 12 * i;
      for(int j = 3; j < 12; j++) {
        ptr[j] = scaleFactor*(double(*((float*)(triRaw + 4 * j))));
      }
      double v1[3], v2[3];
      for(unsigned iDim = 0; iDim < 3; iDim++) {
        v1[iDim] = ptr[iDim + 3 * 2] - ptr[iDim + 3 * 1];
        v2[iDim] = ptr[iDim + 3 * 3] - ptr[iDim + 3 * 2];
      }
      Geo3d::crossproduct(ptr, v1, v2);
      Geo3d::normalize(ptr);
    }
    return true;
  }
  else
  {
    std::cerr << "File " << fname << " openning error!" << std::endl;
    return false;
  }
}

void StlGeo::GetTriIdList(
  std::set<size_t> &triList, 
  const size_t &cellId
) const {
  triList.clear();
  int loc[3];
  gridIndex.CalCellIJK(loc, cellId);
  Range<3> rg;
  rg.SetStart(loc); rg.SetEnd(loc);
  rg.AddBuffer(1);
  for(int carry = rg.IniCarry(loc); carry == 0; carry = rg.AddOne(loc))
  {
    size_t cellNId = gridIndex.CalCellID(loc);
    if(cellNId != -1) {
      for(auto range = cellToTri.equal_range(cellNId); range.first != range.second; range.first++) {
        size_t triId = range.first->second;
        triList.emplace(triId);
      }
    }
  }
}

bool StlGeo::GenBoundingBox(const double &minL) {
  // boundingBox(6); //	first 3 record the bottom left, last 3 record the top right
  const size_t triNum = get_triangleNum();
  if(triNum < 1) {
    return false;
  }
  cellDx = minL;
  double* ptr = get_triangle(0);
  for(int i = 0; i < 3; i++) {
    boundingBox[i] = boundingBox[i + 3] = ptr[i + 3];
  }
  for(int i = 0; i < triNum; i++) {
    for(int j = 1; j < 4; j++) {
      for(int iDim = 0; iDim < 3; iDim++) {
        if(ptr[iDim + 3 * j] < boundingBox[iDim]) boundingBox[iDim] = ptr[iDim + 3 * j];
        if(ptr[iDim + 3 * j] > boundingBox[iDim + 3]) boundingBox[iDim + 3] = ptr[iDim + 3 * j];
      }
    }
    ptr += 12;
  }
  //////////////////////////////////
  int cellNum[3], nodeNum[3];
  for(int iDim = 0; iDim < 3; iDim++) {
    cellNum[iDim] = (int)floor((boundingBox[iDim + 3]-boundingBox[iDim]) / minL) + 5;
    nodeNum[iDim] = cellNum[iDim] + 1;
  }
  for(int iDim = 0; iDim < 3; iDim++) {
    boundingBox[iDim] = boundingBox[iDim] - 2.5 * minL;
    boundingBox[iDim + 3] = boundingBox[iDim] + cellNum[iDim] * minL;
  }
  nodeRange = Range<3>(0,0,0, nodeNum[0], nodeNum[1], nodeNum[2]);
  gridIndex = RangeIndex<3>(nodeRange, RangeIndex<3>::inputNodeRange);
  std::cout << "StlGeo::fileName : " << fileName << std::endl;
  std::cout << "StlGeo::boundingBox : \n\t" 
    << boundingBox[0] << ' ' << boundingBox[1] << ' ' << boundingBox[2] << "\n\t"
    << boundingBox[3] << ' ' << boundingBox[4] << ' ' << boundingBox[5] << std::endl;
  std::cout << "StlGeo::cellNum : \n\t" 
    << cellNum[0] << ' ' << cellNum[1] << ' ' << cellNum[2] << std::endl;

  int loc[3];
  //size_t cellId;
  // triangle to cell
  std::set<std::pair<size_t,size_t>> triToCellTemp;
  std::set<std::pair<size_t,size_t>>::iterator iter = triToCellTemp.begin();
  std::list<std::array<double,10>> triQue;
  for(size_t triId = 0; triId < triNum; triId++) {
    const double *triNorm = get_triangle(triId);
    const double *triPos = triNorm + 3;
    double maxL = 0.0;
    double tempV[3];

    for(int iEdge = 0; iEdge < 3; iEdge++) {
      int iNode1 = (iEdge + 1) % 3;
      for(int iDim = 0; iDim < 3; iDim++) {
        tempV[iDim] = triPos[iDim + 3 * iEdge] - triPos[iDim + 3 * iNode1];
      }
      double length = sqrt(Geo3d::dotproductVV(tempV,tempV));
      if(length > maxL) maxL = length;
    }

    std::array<double, 10> triElem;
    memcpy(triElem.data(), triPos, 9 * sizeof(double));
    triElem[9] = maxL;
    triQue.push_front(triElem);
    while(triQue.size() != 0) {
      std::array<double, 10> subTri = triQue.front();
      triQue.pop_front();
      if(subTri[9] > 0.99 * minL) {
        double fourTri[3 * 3 * 4];	 
        Geo3d::SplitTriangle(fourTri, subTri.data());
        for(int iTri = 0; iTri < 4; iTri++) {
          std::array<double, 10> tempTri;
          memcpy(tempTri.data(), fourTri + 9 * iTri, 9 * sizeof(double));
          tempTri[9] = subTri[9] * 0.5;
          triQue.push_front(tempTri);
        }
      }
      else {
        Box bb = Box::FindBoundingBoxOf(subTri.data(), 3, 3);
        //std::cout << bb[0] << ' ' << bb[1] << ' ' << bb[2] << "   " << bb[3] << ' ' << bb[4] << ' ' << bb[5] << std::endl;
        std::set<size_t> cellIdSet;
        int rg[6];
        for(int iDim = 0; iDim < 3; iDim++) {
          rg[iDim + 0] = FindCellIndex(bb[iDim + 0], boundingBox[iDim], minL);
          rg[iDim + 3] = FindCellIndex(bb[iDim + 3], boundingBox[iDim], minL);
        }
        //std::cout << rg[0] << ' ' << rg[1] << ' ' << rg[2] << "   " << rg[3] << ' ' << rg[4] << ' ' << rg[5] << std::endl;
        for(int i = 0; i < 2; i++) 
        for(int j = 0; j < 2; j++) 
        for(int k = 0; k < 2; k++) {
          loc[0] = rg[0 + 3 * i];
          loc[1] = rg[1 + 3 * j];
          loc[2] = rg[2 + 3 * k];
          size_t cellId = gridIndex.CalCellID(loc);
          cellIdSet.insert(cellId);
        }
        for(auto cellId : cellIdSet) {
          double localCell[6];
          gridIndex.CalCellIJK(loc, cellId);
          for(int iDim = 0; iDim < 3; iDim++) {
            localCell[iDim] = boundingBox[iDim] + loc[iDim] * minL;
            localCell[iDim + 3] = boundingBox[iDim] + (loc[iDim] + 1) * minL;
          }
          if(Geo3d::is_intersected(subTri.data(), localCell, 1e-10)) 
            iter = triToCellTemp.insert(iter, std::pair<size_t,size_t>(triId, cellId));
          else continue;
        }
      }
    }
  }

  Range<3> crg = gridIndex.GetCellRange();
  Range<3> nrg = gridIndex.GetNodeRange();
  // cell to triangle; triangle to cell
  colorMap = MaskRange<int, 3>(crg, int(geo::ndef));

  for(auto elem : triToCellTemp) {
    //const int triId = elem.first;
    //const int nodeId = elem.second;
    triToCell.emplace(elem.first, elem.second);
    cellToTri.emplace(elem.second, elem.first);
    colorMap.SetMask(elem.second, nnwlbm::geo::boundary);
  }

  //! color map
  while(GenConnectedDomain()) {}

  //WriteColorMap();
  return true;
}

bool StlGeo::GenConnectedDomain() {
  auto crg = gridIndex.GetCellRange();
  int loc[3];
  size_t cellId;
  bool found = false;
  // find first ndef in the colorMap
  for(int carry = crg.IniCarry(loc); carry == 0; carry = crg.AddOne(loc))
  {
    cellId = gridIndex.CalCellID(loc);
    if(colorMap.GetMask(cellId) == geo::ndef) {
      found = true;
      break;
    }
  }
  if(!found) {
    std::cout << "StlGeo::_connected_domain() : all connected domains are found out.\n"
      << "numOfConnectedDomain = "<< numOfConnectedDomain << std::endl;
    return false;
  }

  numOfConnectedDomain++;
  colorMap.SetMask(cellId, numOfConnectedDomain);
  std::cout << "StlGeo::_connected_domain() " << numOfConnectedDomain << std::endl;
  std::queue<size_t, std::list<size_t>> searchQue;

  searchQue.push(cellId);

  double err = 1e-10;
  size_t neighborId;
  int neighborLoc[3];
  while(searchQue.size() != 0) {
    cellId = searchQue.front();
    gridIndex.CalCellIJK(loc, cellId);
    // search neighbor
    for(int iDim = 0; iDim < 3; iDim++) 
      for(int dd = -1; dd < 2; dd += 2) {
        memcpy(neighborLoc, loc, 3*sizeof(int));
        neighborLoc[iDim] += dd;
        if(!crg.IsContaining(neighborLoc)) continue;
        neighborId = gridIndex.CalCellID(neighborLoc);
        if(colorMap.GetMask(neighborId) != nnwlbm::geo::ndef) continue;
        searchQue.push(neighborId);
        colorMap.SetMask(neighborId, numOfConnectedDomain);
      }	// end search neighbor
    searchQue.pop();
  }// while(searchQue.size() != 0)
  return found;
}

void StlGeo::WriteColorMap() const {
  auto crg = gridIndex.GetCellRange();
  std::ofstream outFile(fileName+"_colorMap.dat", std::ios::out);
  outFile << "TITLE = The colorMap of " << fileName << "\n";
  outFile << "VARIABLES = X, Y, Z, color\n";
  outFile << "ZONE I = " << crg.length(0) << ", J = " << crg.length(1) 
    << ", K = " << crg.length(2) << "\n" << std::endl;

  int loc[3];
  for(int carry = crg.IniCarry(loc); carry == 0; carry = crg.AddOne(loc)) {
    size_t cellId = gridIndex.CalCellID(loc);
    outFile << boundingBox[0] + (0.5 + loc[0]) * cellDx << ' ' << boundingBox[1] + (0.5 + loc[1]) * cellDx 
      << ' ' << boundingBox[2] + (0.5 + loc[2]) * cellDx << ' ' << colorMap.GetMask(cellId) << ' ' << std::endl;
  }
  std::cout << "StlGeoDT::write " << fileName+"_colorMap.dat" << std::endl;
}


int StlGeo::CalRelativeLocation(const double *pos, const double &tolerance) const {
  int ans = CalRelativeLocationOrigin(pos, tolerance);
  if(ans > 1) return geo::inside;
  else if(ans == 1) return geo::outside;
  return ans;
}

int StlGeo::CalRelativeLocationOrigin(const double *pos, const double &tolerance) const {
  int loc[3];
  double nPos[3];
  double originC[3];
  for(int iDim = 0; iDim < 3; iDim++) {
    loc[iDim] = FindCellIndex(pos[iDim], boundingBox[iDim + 3 * 0], cellDx);
    originC[iDim] = boundingBox[iDim] + 0.5 * cellDx;
  }
  UniformGrid<3> gridC(cellDx, originC);
  Range<3> crg = gridIndex.GetCellRange();
  if(!crg.IsContaining(loc)) return 1; // the first connected domain
  const size_t cellId = gridIndex.CalCellID(loc);

  if(colorMap.GetMask(cellId) != nnwlbm::geo::boundary) {
    return colorMap.GetMask(cellId);
  }

  auto triangles = cellToTri.equal_range(cellId);
  for(;triangles.first != triangles.second; triangles.first++) {
    if (Geo3d::on_triangle(pos, GetTrianglePos(triangles.first->second), tolerance)) {
      return geo::boundary;
    }
  }

  size_t cellNId;
  Range<3> neighborRange;
  double intersection[3];
  std::set<size_t> connectedCell;
  std::list<size_t> toBeExplored;

  // start exploring

  // try to find connection
  neighborRange.SetStart(loc);
  neighborRange.SetEnd(loc);
  neighborRange.AddBuffer(1);
  neighborRange = Range<3>::CalIntersection(crg, neighborRange);
  bool isConnected = true;
  for(int carry = neighborRange.IniCarry(loc); 
    carry == 0; carry = neighborRange.AddOne(loc))
  {
    cellNId = gridIndex.CalCellID(loc);
    gridC.GetPos(nPos, loc);
    isConnected = true;
    auto triangles = cellToTri.equal_range(cellId);
    for(;triangles.first != triangles.second; triangles.first++) {
      int crossRelation = Geo3d::TriangleIntersection(intersection, 
        GetTrianglePos(triangles.first->second), pos, nPos, tolerance);
      if(crossRelation == geo::outside) continue;
      else {
        isConnected = false;
        break;
      }
    }
    if(isConnected) {
      triangles = cellToTri.equal_range(cellNId);
      for(;triangles.first != triangles.second; triangles.first++) {
        int crossRelation = Geo3d::TriangleIntersection(intersection, 
          GetTrianglePos(triangles.first->second), pos, nPos, tolerance);
        if(crossRelation == geo::outside) continue;
        else {
          isConnected = false;
          break;
        }
      }
      if(isConnected) {
        // if connection with fluid or solid is found, return
        if(colorMap.GetMask(cellNId) != geo::boundary) {
          return colorMap.GetMask(cellNId);
        }
        // extend searching area on boundary
        connectedCell.emplace(cellNId);
        if(cellId != cellNId) toBeExplored.push_back(cellNId);
      }
    }
  }

  while(toBeExplored.size() != 0) {
    size_t cellCId = toBeExplored.front();
    toBeExplored.pop_front();
    double cPos[3];
    gridIndex.CalCellIJK(loc, cellCId);
    gridC.GetPos(cPos, loc);
    // search neighbor
    neighborRange.SetStart(loc);
    neighborRange.SetEnd(loc);
    neighborRange.AddBuffer(1);
    neighborRange = Range<3>::CalIntersection(crg, neighborRange);
    for(int carry = neighborRange.IniCarry(loc); 
      carry == 0; carry = neighborRange.AddOne(loc))
    {
      cellNId = gridIndex.CalCellID(loc);
      if(connectedCell.find(cellNId) != connectedCell.end()) continue;
      gridC.GetPos(nPos, loc);
      isConnected = true;
      auto triangles = cellToTri.equal_range(cellCId);
      for(;triangles.first != triangles.second; triangles.first++) {
        int crossRelation = Geo3d::TriangleIntersection(intersection, 
          GetTrianglePos(triangles.first->second), cPos, nPos, tolerance);
        if(crossRelation == geo::outside) continue;
        else {
          isConnected = false;
          break;
        }
      }
      if(isConnected) {
        triangles = cellToTri.equal_range(cellNId);
        for(;triangles.first != triangles.second; triangles.first++) {
          int crossRelation = Geo3d::TriangleIntersection(intersection, 
            GetTrianglePos(triangles.first->second), cPos, nPos, tolerance);
          if(crossRelation == geo::outside) continue;
          else {
            isConnected = false;
            break;
          }
        }
        if(isConnected) {
          // if connection with fluid or solid is found, return
          if(colorMap.GetMask(cellNId) != geo::boundary) {
            return colorMap.GetMask(cellNId);
          }
          connectedCell.emplace(cellNId);
          //if(cellId != cellNId) toBeExplored.push_back(cellNId);
          toBeExplored.push_back(cellNId);
        } //end if connected
      } // end if connected
    }// end traverse neighbors
  } // end while
  
  return nnwlbm::geo::boundary;
}

std::shared_ptr<MaskRange<int, 3>> StlGeo::sharpen_boundary() const {
  auto ans = std::make_shared<MaskRange<int, 3>>(colorMap);
  auto crg = gridIndex.GetCellRange();
  int loc[3];
  double pos[3];

  for(int carry = crg.IniCarry(loc); carry == 0; carry = crg.AddOne(loc))
  {
    size_t cellId = gridIndex.CalCellID(loc);
    if(colorMap.GetMask(cellId) == geo::boundary) {
      pos[0] = boundingBox[0] + (0.5 + loc[0]) * cellDx;
      pos[1] = boundingBox[1] + (0.5 + loc[1]) * cellDx;
      pos[2] = boundingBox[2] + (0.5 + loc[2]) * cellDx;
      //std::cout << loc[0] <<' ' << loc[1] << ' ' << loc[2] << std::endl;
      //if( loc[0] == 26 && loc[1] == 3 && loc[2] == 9) 
      //{
      //  std::cout << loc[0] <<' ' << loc[1] << ' ' << loc[2] << std::endl;
      //}
      ans->SetMask(cellId, CalRelativeLocationOrigin(pos));
    }
  }

  std::ofstream outFile(fileName+"_sharpen_boundary.dat", std::ios::out);
  outFile << "TITLE = The colorMap of " << fileName << "\n";
  outFile << "VARIABLES = X, Y, Z, color\n";
  outFile << "ZONE I = " << crg.length(0) << ", J = " << crg.length(1) 
    << ", K = " << crg.length(2) << "\n" << std::endl;

  for(int carry = crg.IniCarry(loc); carry == 0; carry = crg.AddOne(loc)) {
    size_t cellId = gridIndex.CalCellID(loc);
    outFile << boundingBox[0] + (0.5 + loc[0]) * cellDx << ' ' << boundingBox[1] + (0.5 + loc[1]) * cellDx 
      << ' ' << boundingBox[2] + (0.5 + loc[2]) * cellDx << ' ' << ans->GetMask(cellId) << ' ' << std::endl;
  }
  std::cout << "StlGeoDT::write " << fileName+"_sharpen_boundary.dat" << std::endl;

  return ans;
}

bool StlGeo::cal_normal(double *norm, const double *pos, const double &tolerance) const
{
  return true;
}

bool StlGeo::CalIntersection(double *intersection, const double *p1, const double *p2, 
  const double &tolerance) const
{
  int loc1[3], loc2[3], mLoc[3];
  double unitVec[3], middleP[3];
  for(int iDim = 0; iDim < 3; iDim++) {
    loc1[iDim] = FindCellIndex(p1[iDim], boundingBox[iDim + 3 * 0], cellDx);
    loc2[iDim] = FindCellIndex(p2[iDim], boundingBox[iDim + 3 * 0], cellDx);
    unitVec[iDim] = p2[iDim] - p1[iDim];
  }
  size_t cell1Id = gridIndex.CalCellID(loc1);
  size_t cell2Id = gridIndex.CalCellID(loc2);
  double totalL = Geo3d::normalize(unitVec);
  size_t middleId = cell2Id;
  std::set<size_t> triList;
  GetTriIdList(triList, middleId);
  if(triList.size() == 0) {
    middleId = cell1Id;
    GetTriIdList(triList, middleId);
  }
  for(double dl = 0.0; dl < totalL; ) {
    for(auto triId : triList) {
      const double* triPos = GetTrianglePos(triId);
      const double* triNorm = GetTriangleNorm(triId);
      if(Geo3d::TriangleIntersection(intersection, triPos, p1, p2, tolerance) != nnwlbm::geo::outside) {
        //if(!Geo3d::on_triangle(intersection, triPos)) {
        //	cout << "WTF" << endl;
        //	Geo3d::TriangleIntersection(intersection, triPos, p1, p2);
        //	Geo3d::on_triangle(intersection, triPos);
        //}
        return true;
      }
    }
    dl += 2 * cellDx;
    for(int i = 0; i < 3; i++) {
      middleP[i] = p1[i] + dl * unitVec[i];
      mLoc[i] = FindCellIndex(middleP[i], boundingBox[i + 3 * 0], cellDx);
    }
    middleId = gridIndex.CalCellID(mLoc);
    GetTriIdList(triList, middleId);
  }
  if(triList.size() == 0) return false;
  return false;
}
}
