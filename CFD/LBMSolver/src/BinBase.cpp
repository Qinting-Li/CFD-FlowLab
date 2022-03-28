////////////////////////////////////////////////////////
// author: Chen, Songze. email: jacksongze@hust.edu.cn
////////////////////////////////////////////////////////

#include "CommonLibDll.h"
#include "BinBase.h"

#ifdef _WIN32
#include <io.h>
#include <direct.h> 
#else
#include <unistd.h>
#include <sys/stat.h>
#endif
#include <vector>

#ifdef _WIN32
#define ACCESS(fileName,accessMode) _access(fileName,accessMode)
#define MKDIR(path) _mkdir(path)
#else
#define ACCESS(fileName,accessMode) access(fileName,accessMode)
#define MKDIR(path) mkdir(path,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif

int create_directory(const std::string &directoryPath)
{
  size_t dirPathLen = directoryPath.length();
  std::vector<char> tmpDirPath(dirPathLen+1,'\0');
  for (int i = 0; i < dirPathLen; ++i)
  {
    tmpDirPath[i] = directoryPath[i];
    if (tmpDirPath[i] == '\\' || tmpDirPath[i] == '/')
    {
      if (ACCESS(tmpDirPath.data(), 0) != 0)
      {
        tmpDirPath[i + 1] = '\0';
        int ret = MKDIR(tmpDirPath.data());
        if (ret != 0)
        {
          return ret;
        }
      }
    }
  }
  return 0;
}

bool nnwlbm::DataClass::output_bin(
  const DataClass &ref, /**< [in]需要输出的数据 */
  const std::string &filename,  /**< [in]输出的文件名 */
  const std::string &path /**< [in]输出的路径，默认为当前目录 */
)
{ 
  bool flag = true;
  if(create_directory(path) < 0) return false;
  std::ofstream ostr;
  ostr.open(path+filename, std::ios::out | std::ios::binary);
  if(!ref.OutputHead(ostr, ref.name())) flag = false;
  if(!ref.OutputBinPrivate(ostr,filename,path)) flag = false;
  ostr.close();
  return flag;
}

bool nnwlbm::DataClass::OutputHead(
  std::ofstream &ostr, /**< [in,out]二进制文件输出流 */
  const std::string &s  /**< [in] 文件头的内容，一般为输出类型的名称 */
) 
{
  //ostr.setf( std::ios::binary );
  size_t num = s.size();
  ostr.write((const char*)(&num), sizeof(int));
  ostr.write(s.c_str(), num);
  return true;
}

bool nnwlbm::DataClass::OutputBinPrivate(
  std::ofstream &ostr,  /**< [in,out]二进制文件输出流 */
  const std::string &filename, /**< [in]输出的文件名 */
  const std::string & path /**< [in]输出的路径，默认为当前目录 */
) const
{ return false; }

::std::shared_ptr<nnwlbm::DataClass> nnwlbm::BinReader::read(
  const std::string &filename, /**< [in]文件名 */
  const std::string &path /**< [in]文件路径 */
) const 
{
  std::ifstream istr;
  istr.open(path+filename, std::ios::in | std::ios::binary);
  if(!istr.is_open()) {
    std::cout << "cannot open " << path << filename << std::endl;
    return nullptr;
  }
  std::string type = InputHead(istr);
  inptr ptr = get(type);
  if(ptr == nullptr) {
    istr.close();
    return nullptr;
  }
  else { 
    auto ans = ptr->InputBinPrivate(istr, filename, path);
    istr.close();
    return ans;
  }
}

namespace nnwlbm {
  nnwlbm::BinReader binReader;
}


