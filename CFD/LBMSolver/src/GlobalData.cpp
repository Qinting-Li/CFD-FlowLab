///////////////////////////////////////////////////////////////
// author : He, Jiaming. email: 2457837217@qq.com
///////////////////////////////////////////////////////////////
// participant 1 : Chen, Songze. email: jacksongze@hust.edu.cn
///////////////////////////////////////////////////////////////

//GlobalData.cpp
#include "CommonLibDll.h"
#include "GlobalData.h"
#include <vector>

namespace nnwlbm {

bool GlobalData::OutputBinPrivate(
  std::ofstream &ostr,
  const std::string &filename,
  const std::string &path
)const
{
  std::string memberPath = path + filename + "_g\\";
  size_t num = this->name().size();
  ostr.write((const char*)(&num), sizeof(size_t));
  ostr.write(this->name().c_str(), num);
  for (auto elem : this->dataMap) {
    if (elem.second == nullptr) continue;
    //if (std::dynamic_pointer_cast<DataClass>(elem.second) != nullptr) continue;
    size_t num = elem.first.size();
    ostr.write((const char*)(&num), sizeof(size_t));
    ostr.write(elem.first.c_str(), elem.first.size() * sizeof(char));
    elem.second->output_bin(*(elem.second), elem.first, memberPath);
  }
  return true;
}

::std::string GlobalData::name() const { return "GlobalData"; }

template<>
Dptr BinIn<GlobalData>::InputBinPrivate(
  std::ifstream &istr,
  const std::string &filename,
  const std::string &path
) const
{
  std::vector<char> tempMem(sizeof(int));
  istr.read(tempMem.data(), sizeof(int));
  unsigned nameLength = *((unsigned*)tempMem.data());
  tempMem.resize(nameLength);
  istr.read(tempMem.data(), nameLength);
  auto ans = std::make_shared<GlobalData>();
  if (ans == nullptr) return ans;
  std::string memberPath = path + filename + "_g\\";
  for (istr.read(tempMem.data(), sizeof(int)); !istr.eof(); istr.read(tempMem.data(), sizeof(int)))
  {
    int tagLength = *((int*)tempMem.data());
    tempMem.resize(tagLength);
    istr.read(tempMem.data(), tagLength);
    std::string tag(tempMem.data(), tagLength);
    Dptr member = binReader.read(tag, memberPath);
    bool isbool;
    isbool = ans->emplace(tag, member);
  }
  return ans;
}

template<>
::std::string BinIn<GlobalData>::name() const { return "GlobalData"; }

}

// end of GlobalData.cpp
