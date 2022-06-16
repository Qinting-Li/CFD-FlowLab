
#include "CommonLibDll.h"
#include "PHLBMUnstructGrid.h"

namespace phlbm {

  //std::vector<double> gen_triangleVec_from() {
  std::vector<double> gen_triangleVec_from_region() {
    PHSPACE::Region region;
    region.ReadGrid();
    PHSPACE::Zone* zPtr = region.GetZone(0);
    PHSPACE::Grid* phGrid = zPtr->GetGeometry()->GetGrid();
    PHSPACE::UnstructGrid *unsPhGrid = PHSPACE::UnstructGridCast(phGrid);

    int nBoundFace = unsPhGrid->GetNBoundFace();
    int *face2node = unsPhGrid->GetFace2Node();
    int *nodeNumberOfEachFace  = unsPhGrid->GetNodeNumberOfEachFace();

    PHSPACE::RDouble *x = unsPhGrid->GetX();
    PHSPACE::RDouble *y = unsPhGrid->GetY();
    PHSPACE::RDouble *z = unsPhGrid->GetZ();

    //! Boundary faces.
    int nodePosition = 0;
    int triNum = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
      triNum += nodeNumberOfEachFace[iFace] - 2;
      nodePosition += nodeNumberOfEachFace[iFace];
    }
    std::vector<double> ans(12 * triNum);
    nodePosition = 0;
    triNum = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
      for (int jNode = 2; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
      {
        double* tri = ans.data() + 12 * triNum;
        int id = face2node[nodePosition];
        tri[0 + 3 * 1] = x[id];
        tri[1 + 3 * 1] = y[id];
        tri[2 + 3 * 1] = z[id];

        id = face2node[nodePosition + jNode - 1];
        tri[0 + 3 * 2] = x[id];
        tri[1 + 3 * 2] = y[id];
        tri[2 + 3 * 2] = z[id];
        id = face2node[nodePosition + jNode];
        tri[0 + 3 * 3] = x[id];
        tri[1 + 3 * 3] = y[id];
        tri[2 + 3 * 3] = z[id];

        double v1[3], v2[3];
        for(unsigned iDim = 0; iDim < 3; iDim++) {
          v1[iDim] = tri[iDim + 3 * 2] - tri[iDim + 3 * 1];
          v2[iDim] = tri[iDim + 3 * 3] - tri[iDim + 3 * 2];
        }
        Geo3d::crossproduct(tri, v1, v2);
        Geo3d::normalize(tri);
        triNum++;
      }
      nodePosition += nodeNumberOfEachFace[iFace];
    }
    return std::move(ans);
  }
}

