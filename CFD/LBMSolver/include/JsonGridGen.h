#pragma once

#include "JsonReader.h"
#include "UniformGrid.h"

namespace phlbm 
{
  template<size_t dim>
  std::shared_ptr<nnwlbm::DataClass> nnwlbm::Json(
    const phlbm::NamedDataBlock &ndb
  ) {
    if(ndb.get_type() != tag::grid) return nullptr;
    int dimVar;
    if(!ndb.GetData("dim", &dimVar, 1)) {
      std::ostringstream oss;
      oss << "nnwlbm::UniformGrid<dim> needs to specify \"dim\" key in DataBlock.";
      PHSPACE::TK_Exit::ExceptionExit(oss);
    }
    if(dimVar != dim) {
      PHSPACE::TK_Exit::UnexpectedVarValue("dimVar in gen_grid_t", dimVar);
    }
    double origin[dim];
    double dx;
    if(ndb.GetData("origin", origin, dim) 
      && ndb.GetData("dx", &dx, 1)) {
      return std::make_shared<nnwlbm::UniformGrid<dim>>(dx, origin);
    }
    else return nullptr;
  }
}
