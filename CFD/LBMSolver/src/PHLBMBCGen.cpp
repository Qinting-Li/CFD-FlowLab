// PHLBMBCGen.cpp
#include <CommonLibDll.h>
#include <PHLBMGenerator.h>
#include <Geo_SimpleBC.h>
#include <BoundaryFactory.h>

namespace phlbm {
template<>
Generator<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>>::Generator() {
  EnrollInput<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>,
    nnwlbm::BounceBackFactory<nnwlbm::D2Q9>,
    nnwlbm::MovingBounceBackFactory<nnwlbm::D2Q9>,
    nnwlbm::NonEquilibriumExtrapolationFactory<nnwlbm::D2Q9>
  >(inputMap);
}

template<>
::std::shared_ptr<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>> Generator<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>>::operator()(
  const NamedDataBlock& ndb, // [in,out]生成几何的数据
  const std::map<std::string, NamedDataBlock>& ndbMap // [in,out]所有经风雷导入的数据 
  ) const {
  ndb.get_name();
  if (ndb.get_type() != tag::bc && ndb.get_type() != -tag::bc) return nullptr;
  string bcTypeName;
  ndb.GetData("bcTypeName", &bcTypeName, 1);
  auto generator = inputMap.get(bcTypeName);
  if (generator == nullptr) return nullptr;
  return generator->InputPrivate(ndb, ndbMap);
}

template<>
::std::string In<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>, nnwlbm::BounceBackFactory<nnwlbm::D2Q9>>::name() const { return "BounceBackFactory<D2Q9>"; }

template<>
std::shared_ptr<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>> // return value
In<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>, nnwlbm::BounceBackFactory<nnwlbm::D2Q9>>::InputPrivate(
  const NamedDataBlock& ndb,
  const std::map<std::string, NamedDataBlock>& ndbMap
) const {
  return std::make_shared<nnwlbm::BounceBackFactory<nnwlbm::D2Q9>>();
}

template<>
::std::string In<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>, nnwlbm::MovingBounceBackFactory<nnwlbm::D2Q9>>::name() const { return "MovingBounceBackFactory<D2Q9>"; }

template<>
std::shared_ptr<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>> // return value
In<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>, nnwlbm::MovingBounceBackFactory<nnwlbm::D2Q9>>::InputPrivate(
  const NamedDataBlock& ndb,
  const std::map<std::string, NamedDataBlock>& ndbMap
) const {
  double latticeVel;
  ndb.GetData("latticeVel", &latticeVel, 1);
  double macroVar[2];
  ndb.GetData("macroVar", macroVar, 2);
  return std::make_shared<nnwlbm::MovingBounceBackFactory<nnwlbm::D2Q9>>(latticeVel, macroVar);
}

template<>
::std::string In<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>, nnwlbm::NonEquilibriumExtrapolationFactory<nnwlbm::D2Q9>>::name() const { return "NonEquilibriumExtrapolationFactory<D2Q9>"; }

template<>
std::shared_ptr<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>> // return value
In<nnwlbm::BoundaryFactory<nnwlbm::D2Q9>, nnwlbm::NonEquilibriumExtrapolationFactory<nnwlbm::D2Q9>>::InputPrivate(
  const NamedDataBlock &ndb,
  const std::map<std::string, NamedDataBlock> &ndbMap
) const {
  double latticeVel;
  ndb.GetData("latticeVel", &latticeVel, 1);
  double macroVar[2];
  ndb.GetData("macroVar", macroVar, 2);
  int vectorVar[2];
  ndb.GetData("vector", vectorVar, 2);
  return std::make_shared<nnwlbm::NonEquilibriumExtrapolationFactory<nnwlbm::D2Q9>>(latticeVel, macroVar, vectorVar);
}


template<>
Generator<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>>::Generator() {
  EnrollInput<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>,
    nnwlbm::BounceBackFactory<nnwlbm::D3Q15>,
    nnwlbm::MovingBounceBackFactory<nnwlbm::D3Q15>,
    nnwlbm::NonEquilibriumExtrapolationFactory<nnwlbm::D3Q15>,
    nnwlbm::NonEquilibriumOfPressureFactory<nnwlbm::D3Q15>
  >(inputMap);
}

template<>
::std::shared_ptr<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>> Generator<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>>::operator()(
  const NamedDataBlock &ndb, // [in,out]生成几何的数据
  const std::map<std::string, NamedDataBlock> &ndbMap // [in,out]所有经风雷导入的数据 
  ) const {
  ndb.get_name();
  if (ndb.get_type() != tag::bc && ndb.get_type() != -tag::bc) return nullptr;
  string bcTypeName;
  ndb.GetData("bcTypeName", &bcTypeName, 1);
  auto generator = inputMap.get(bcTypeName);
  if (generator == nullptr) return nullptr;
  return generator->InputPrivate(ndb, ndbMap);
}

template<>
::std::string In<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>, nnwlbm::BounceBackFactory<nnwlbm::D3Q15>>::name() const { return "BounceBackFactory<D3Q15>"; }

template<>
std::shared_ptr<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>> // return value
In<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>, nnwlbm::BounceBackFactory<nnwlbm::D3Q15>>::InputPrivate(
  const NamedDataBlock &ndb,
  const std::map<std::string, NamedDataBlock> &ndbMap
) const {
  return std::make_shared<nnwlbm::BounceBackFactory<nnwlbm::D3Q15>>();
}

template<>
::std::string In<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>, nnwlbm::MovingBounceBackFactory<nnwlbm::D3Q15>>::name() const { return "MovingBounceBackFactory<D3Q15>"; }

template<>
std::shared_ptr<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>> // return value
In<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>, nnwlbm::MovingBounceBackFactory<nnwlbm::D3Q15>>::InputPrivate(
  const NamedDataBlock &ndb,
  const std::map<std::string, NamedDataBlock> &ndbMap
) const {
  double latticeVel;
  ndb.GetData("latticeVel", &latticeVel, 1);
  double macroVar[3];
  ndb.GetData("macroVar", macroVar, 3);
  return std::make_shared<nnwlbm::MovingBounceBackFactory<nnwlbm::D3Q15>>(latticeVel, macroVar);
}

template<>
::std::string In<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>, nnwlbm::NonEquilibriumExtrapolationFactory<nnwlbm::D3Q15>>::name() const { return "NonEquilibriumExtrapolationFactory<D3Q15>"; }

template<>
std::shared_ptr<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>> // return value
In<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>, nnwlbm::NonEquilibriumExtrapolationFactory<nnwlbm::D3Q15>>::InputPrivate(
  const NamedDataBlock &ndb,
  const std::map<std::string, NamedDataBlock> &ndbMap
) const {
  double latticeVel;
  ndb.GetData("latticeVel", &latticeVel, 1);
  double macroVar[3];
  ndb.GetData("macroVar", macroVar, 3);
  int vectorVar[3];
  ndb.GetData("vector", vectorVar, 3);
  return std::make_shared<nnwlbm::NonEquilibriumExtrapolationFactory<nnwlbm::D3Q15>>(latticeVel, macroVar, vectorVar);
}

template<>
::std::string In<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>, nnwlbm::NonEquilibriumOfPressureFactory<nnwlbm::D3Q15>>::name() const { return "NonEquilibriumOfPressureFactory<D3Q15>"; }

template<>
std::shared_ptr<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>> // return value
In<nnwlbm::BoundaryFactory<nnwlbm::D3Q15>, nnwlbm::NonEquilibriumOfPressureFactory<nnwlbm::D3Q15>>::InputPrivate(
  const NamedDataBlock &ndb,
  const std::map<std::string, NamedDataBlock> &ndbMap
) const {
  double latticeVel;
  ndb.GetData("latticeVel", &latticeVel, 1);
  double pressure;
  ndb.GetData("pressure", &pressure, 1);
  int vectorVar[3];
  ndb.GetData("vector", vectorVar, 3);
  return std::make_shared<nnwlbm::NonEquilibriumOfPressureFactory<nnwlbm::D3Q15>>(latticeVel, pressure, vectorVar);
}

}
//
//////////////////////////////////////////////////
////
////    请添加各种边界条件生成工具
////
//////////////////////////////////////////////////
