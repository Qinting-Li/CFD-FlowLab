///////////////////////////////////////////////////////////////
// author : He, Jiaming. email: 2457837217@qq.com
///////////////////////////////////////////////////////////////
// participant 1 : Wang,Junquan. email: as2697920806@gmail.com
///////////////////////////////////////////////////////////////
//GetConfig.cpp
#include "CommonLibDll.h"
#include "GetConfig.h"
bool nnwlbm::GetConfig::is_space(
  const char &c
)
{
  if (c == ' ' || c == '\t')
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool nnwlbm::GetConfig::is_commentChar(
  const char &c
)const
{
  if (c == '#')
  {
    return true;
  }
  else
  {
    return false;
  }
}

void nnwlbm::GetConfig::trim(
  std::string &str
)const
{
  if (str.empty())
  {
    return;
  }
  size_t i;
  size_t startPos;
  size_t endPos;
  for (i = 0; i < str.size(); i++)
  {
    if (!nnwlbm::GetConfig::is_space(str[i]))
    {
      break;
    }
  }
  if (i == str.size())
  {
    str = "";
    return;
  }
  startPos = i;

  for (i = str.size() - 1; i >= 0; i--)
  {
    if (!nnwlbm::GetConfig::is_space(str[i]))
    {
      break;
    }
    if(i == 0) break;
  }
  endPos = i;
  str = str.substr(startPos, endPos - startPos + 1);

}

bool nnwlbm::GetConfig::AnalyseLine(
  const std::string &line,
  std::string &key,
  std::string &value
)const
{
  if (line.empty())
  {
    return false;
  }
  size_t startPos = 0, endPos = line.size() - 1, pos;
  if ((pos = line.find('#')) != std::string::npos)
  {
    if (0 == pos)
    {
      return false;
    }
    endPos = pos - 1;
  }
  std::string new_line = line.substr(startPos, endPos - startPos + 1);
  if ((pos = new_line.find("=")) == std::string::npos)
  {
    return false;
  }
  key = new_line.substr(0, pos);
  value = new_line.substr(pos + 1, endPos + 1 - (pos + 1));
  trim(key);
  if (key.empty())
  {
    return false;
  }
  trim(value);
  return true;
}

bool nnwlbm::GetConfig::ReadConfig(
  const std::string &filename
)
{
  ConfigMap.clear();
  std::ifstream inFile(filename.c_str());
  if (!inFile)
  {
    std::cout << "file open failed!" << std::endl;
    return false;
  }
  std::string line, key, value;
  while (std::getline(inFile, line))
  {
    if (AnalyseLine(line, key, value))
    {
      double value2 = atof(const_cast<const char*>(value.c_str()));
      ConfigMap[key] = value2;
    }
  }
  inFile.close();
  return true;
}

void nnwlbm::GetConfig::PrintConfig()const
{
  std::map<std::string, double>::const_iterator mite;
  for (mite = ConfigMap.begin(); mite != ConfigMap.end(); mite++)
  {
    std::cout << mite->first << "=" << mite->second << std::endl;
  }
}


double nnwlbm::GetConfig::FindInConfig(
  const std::string &key
)
{
  std::map<std::string, double>::iterator it;
  it = ConfigMap.find(key);
  if (it == ConfigMap.end())
  {
    std::cout << "FindInConfig is fail";
    return -10000000;
  }
  else
  {
    return it->second;
  }
}
// end of GetConfig.cpp