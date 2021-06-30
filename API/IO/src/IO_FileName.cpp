#include "IO_FileName.h"
#include "TK_Exit.h"
#include "TK_Parse.h"

namespace PHSPACE
{
void GetFileNameofThisProcessor(string &filename)
{
    //! If serial run, do not change file name.
    if (!PHMPI::IsParallelRun())
    {
        return;
    }

    string fname, fext;
    GetNameExt(filename, fname, fext, ".");
	int myID = PHMPI::GetCurrentProcessorID();
    filename = AddDigitToFileName(fname, fext, ".", myID);
}

string AddDigitToFileName(const string &fname, const string &fext, const string &sep, int num)
{
    std::ostringstream oss;
    if (fname.empty()) return oss.str();

    oss << fname << num << sep << fext;
    return oss.str();
}

string GetDefaultDirectory()
{
    return "./results/";
}

string AbsoluteFilename(const string &relativeFilename)
{
    string filename = GetDefaultDirectory();
    filename += relativeFilename;
    return filename;
}

void GetFileNameExtension(const string &fullName, string &mainName, string &extensionName, const string &fileNameSeparator)
{
    basic_string <char>::size_type index;

    index         = fullName.find_last_of(fileNameSeparator);
    mainName      = fullName.substr(0, index);
    extensionName = fullName.substr(index+1, fullName.length() - index - 1);
}

string ChangeExtensionOfFileName(const string &fileName, const string &newExtension)
{
    string mainName, extensionName;
    GetFileNameExtension(fileName, mainName, extensionName, ".");

    ostringstream oss;
    oss << mainName << "." << newExtension;
    return oss.str();
}

int ExtractNumberOfPartFromFileName(const string & source)
{
    int prePosition = -1;
    for (int iChar = 1; iChar < static_cast<int>(source.size()); ++ iChar)
    {
        char ch = source[iChar];
        if (ch == '_')
        {
            if (ch == source[iChar - 1])
            {
                prePosition = iChar;
                break;
            }
        }
    }

    if (prePosition == -1)
    {
        TK_Exit::PrintDebugInfoExit("Error: Can not Extract Number Of Parts From File Name!");
        exit(0);
    }

    string number;
    for (std::size_t iChar = prePosition + 1; iChar < source.size(); ++ iChar)
    {
        char ch = source[iChar];
        if (ch != '_' && ch != '-' && ch != '.')
        {
            number.push_back(ch);
        }
        else
        {
            break;
        }
    }

    int numberParts = -1;
    from_string<int>(numberParts, number, std::dec);
    return numberParts;
}

}