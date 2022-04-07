#include <CommonLibDll.h>
#include <PHLBMGenerator.h>
#include <LBMGeometry.h>
#include "StlGeo.h"

namespace phlbm {
template<>
Generator<nnwlbm::Geometry>::Generator(){
  EnrollInput<nnwlbm::Geometry,
    nnwlbm::EmptyGeo,
    nnwlbm::Box,
    nnwlbm::Circle,
    nnwlbm::UnionGeo,
    nnwlbm::StlGeo
  >(inputMap);
}

template<>
::std::shared_ptr<nnwlbm::Geometry> Generator<nnwlbm::Geometry>::operator()(
  const NamedDataBlock &ndb, // [in,out]生成几何的数据
  const std::map<std::string, NamedDataBlock> &ndbMap // [in,out]所有经风雷导入的数据 
) const {
  ndb.get_name();
  if(ndb.get_type() != tag::geo && ndb.get_type() != -tag::geo) return nullptr;
  string geoTypeName;
  ndb.GetData("geoTypeName", &geoTypeName, 1);
  auto generator = inputMap.get(geoTypeName);
  if(generator == nullptr) {
    std::cout << "Cannot find generator for \"" << geoTypeName << "\"" << std::endl;
    return nullptr;
  }
  return generator->InputPrivate(ndb, ndbMap);
}

template<>
::std::string In<nnwlbm::Geometry, nnwlbm::EmptyGeo>::name() const { return "EmptyGeo"; }
template<>
std::shared_ptr<nnwlbm::Geometry> // return value
In<nnwlbm::Geometry, nnwlbm::EmptyGeo>::InputPrivate(
  const NamedDataBlock &ndb, 
  const std::map<std::string, NamedDataBlock> &ndbMap
) const { 
  int dim = 0;
  ndb.GetData("dim", &dim, 1);
  return std::make_shared<nnwlbm::EmptyGeo>(dim); 
}

template<>
::std::string In<nnwlbm::Geometry, nnwlbm::Box>::name() const { return "Box"; }
template<>
std::shared_ptr<nnwlbm::Geometry> // return value
In<nnwlbm::Geometry, nnwlbm::Box>::InputPrivate(
  const NamedDataBlock &ndb, const std::map<std::string, NamedDataBlock> &ndbMap
) const { 
  int dim = 0;
  ndb.GetData("dim", &dim, 1);
  std::vector<double> box(2 * dim);
  ndb.GetData("box", box.data(), 2 * dim);
  return std::make_shared<nnwlbm::Box>(box); 
}

template<>
::std::string In<nnwlbm::Geometry, nnwlbm::Circle>::name() const { return "Circle"; }
template<>
std::shared_ptr<nnwlbm::Geometry>
In<nnwlbm::Geometry, nnwlbm::Circle>::InputPrivate(
  const NamedDataBlock &ndb, const std::map<std::string, NamedDataBlock> &ndbMap
) const {
  int dim = 0;
  ndb.GetData("dim", &dim, 1);
  double center[2];
  ndb.GetData("center", center, dim);
  double radius;
  ndb.GetData("radius", &radius, 1);
  return std::make_shared<nnwlbm::Circle>(center, radius);
}

template<>
::std::string In<nnwlbm::Geometry, nnwlbm::UnionGeo>::name() const { return "UnionGeo"; }
template<>
std::shared_ptr<nnwlbm::Geometry>
In<nnwlbm::Geometry, nnwlbm::UnionGeo>::InputPrivate(
  const NamedDataBlock &ndb, const std::map<std::string, NamedDataBlock> &ndbMap
) const {
  Generator<nnwlbm::Geometry> geoGen;
  std::string string1;
  ndb.GetData("geometry1", &string1, 1);
  auto geometry1 = geoGen(ndbMap.at(string1), ndbMap);
  std::string string2;
  ndb.GetData("geometry2", &string2, 1);
  auto geometry2 = geoGen(ndbMap.at(string2), ndbMap);
  if (geometry1->GetDim() != geometry2->GetDim()) return nullptr;
  return std::make_shared<nnwlbm::UnionGeo>(geometry1, geometry2);
}

template<>
::std::string In<nnwlbm::Geometry, nnwlbm::StlGeo>::name() const { return "StlGeo"; }
template<>
std::shared_ptr<nnwlbm::Geometry>
In<nnwlbm::Geometry, nnwlbm::StlGeo>::InputPrivate(
  const NamedDataBlock &ndb, const std::map<std::string, NamedDataBlock> &ndbMap
) const {
  auto stlGeo = std::make_shared<nnwlbm::StlGeo>();
  std::string stlName;
  ndb.GetData("stlName", &stlName, 1);
  stlGeo->ReadFile(stlName);
  return stlGeo;
}

}

////////////////////////////////////////////////
//
//    请添加各种几何体的生成工具
//
////////////////////////////////////////////////
