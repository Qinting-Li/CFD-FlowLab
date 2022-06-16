#pragma once


#include "PHLBM.h"
#include "TK_Exit.h"
#include "Range.h"

namespace phlbm 
{
  template<size_t dim>
  std::shared_ptr<nnwlbm::Range<dim>> GenRangeT(
    const phlbm::NamedDataBlock &ndb
  ) {
    if(ndb.get_type() != tag::range) return nullptr;
    int dimVar;
    if(!ndb.GetData("dim", &dimVar, 1)) {
      std::ostringstream oss;
      oss << "phlbm::GenRangeT<dim>(ndb) needs to specify \"dim\" key in DataBlock(ndb).";
      PHSPACE::TK_Exit::ExceptionExit(oss);
    }
    if(dimVar != dim) {
      PHSPACE::TK_Exit::UnexpectedVarValue("phlbm::GenRangeT<dim>(ndb) dim in the files", dimVar);
    }
    int rg[2*dim];
    if(ndb.GetData("range", rg, 2 * dim)) {
      return std::make_shared<nnwlbm::Range<dim>>(rg);
    }
    else return nullptr;
  }
}
