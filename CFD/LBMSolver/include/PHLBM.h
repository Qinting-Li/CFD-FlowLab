#pragma once

#include <mpi.h>
#include <TK_Log.h>
#include <Geo_SimpleBC.h>
#include <Constants.h>
#include <map>

namespace phlbm {

  enum tag {
    geo   = 1,
    bc    = 2,
    grid  = 3,
    range = 4
  };
  class NamedDataBlock {
     PHSPACE::SimpleBC &dblk;
  public:
    NamedDataBlock( PHSPACE::SimpleBC &bc) : dblk(bc){}
    ~NamedDataBlock() {}
    inline int get_type() const { return dblk.GetBCType(); }
    inline const string &get_name() const { return dblk.GetBCName(); }
    inline bool GetData(const string &name, int *data, int size = 1) const
    {
      if(dblk.GetBCParamDataBase()->IsExist(name, PHSPACE::PHINT, size)) { 
        dblk.GetData(name, (void*)data, PHSPACE::PHINT, size); 
        return true;
      }
      PHSPACE::WriteLogFile("Cannot get int data with name \""+name+".\"");
      return false;
    }
    inline bool GetData(const string &name, double *data, int size = 1) const
    { 
      if(dblk.GetBCParamDataBase()->IsExist(name, PHSPACE::PHDOUBLE, size)) { 
        dblk.GetData(name, (void*)data, PHSPACE::PHDOUBLE, size);
        return true; 
      }
      PHSPACE::WriteLogFile("Cannot get double data with name \""+name+".\"");
      return false;
    }
    inline bool GetData(const string &name, float *data, int size = 1) const
    {
      if(dblk.GetBCParamDataBase()->IsExist(name, PHSPACE::PHFLOAT, size)) { 
        dblk.GetData(name, (void*)data, PHSPACE::PHFLOAT, size); 
        return true;
      }
      PHSPACE::WriteLogFile("Cannot get float data with name \""+name+".\"");
      return false;
    }
    inline bool GetData(const string &name, string *data, int size = 1) const
    {
      if(dblk.GetBCParamDataBase()->IsExist(name, PHSPACE::PHSTRING, size)) { 
        dblk.GetData(name, (void*)data, PHSPACE::PHSTRING, size); 
        return true;
      }
      PHSPACE::WriteLogFile("Cannot get string data with name \""+name+".\"");
      return false;
    }
    PHSPACE::Data_Param *GetDataPack() const { return dblk.GetBCParamDataBase(); }
  };

  class BlockTextReader : private PHSPACE::GlobalBoundaryCondition {
    typedef PHSPACE::SimpleBC _phbc;
  public:
    static inline size_t GetNumOfDataBlock() { return GetNumberofBC(); }
    static inline NamedDataBlock GetDataBlock(const int &ib) { return NamedDataBlock(*GetBC(ib)); }
    static inline std::vector<NamedDataBlock> GetDataBlockVector() 
    { 
      std::vector<NamedDataBlock> ans;
      if(GetNumOfDataBlock() == 0) return std::move(ans);
      std::vector<PHSPACE::SimpleBC*>* bcVec = GetGlobalBoundaryConditionList(); 

      for(auto bcPtr : *bcVec) { 
        if(bcPtr != nullptr) ans.push_back(NamedDataBlock(*bcPtr)); 
      }
      return std::move(ans);
    }
    static inline std::map<std::string, NamedDataBlock> GetDataBlockMap() 
    { 
      std::map<std::string, NamedDataBlock> ans;
      if(GetNumOfDataBlock() == 0) return std::move(ans);
      std::vector<PHSPACE::SimpleBC*> *bcVec = GetGlobalBoundaryConditionList(); 

      for(auto bcPtr : *bcVec) { 
        if(bcPtr != nullptr) ans.emplace(bcPtr->GetBCName(), NamedDataBlock(*bcPtr)); 
      }
      return std::move(ans);
    }
    static inline void ReadFromFile() { ReadLBMPara(); }
    //static inline void ReadFromFile(const string& filename) { ReadGlobalBoundaryCondition(); }
  };
}


