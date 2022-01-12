#include "GlobalDataBase.h"
#include "Constants.h"
#include "DataContainer.h"

using namespace std;
#pragma warning (disable:913)
namespace PHSPACE
{
Data_Param * GlobalDataBase::globalParam = new Data_Param();
Data_Field * GlobalDataBase::globalField = new Data_Field();

void GlobalDataBase::UpdateData(const string &name, void *data, int type, int size)
{
    globalParam->UpdateData(name, data, type, size);
}

void GlobalDataBase::GetData(const string &name, void *data, int type, int size)
{
    globalParam->GetData(name, data, type, size);
}

int GlobalDataBase::GetSizeOfData(const string &name)
{
    return globalParam->GetSizeOfData(name);
}

void GlobalDataBase::UpdateDataPtr(const string &name, void *data)
{
    globalField->UpdateDataPtr(name, data);
};

void * GlobalDataBase::GetDataPtr(const string &name)
{
    return globalField->GetDataPtr(name);
};

set < Data_SafeData > * GlobalDataBase::GetBaseP()
{
    return globalParam->GetBaseP();
};

Data_Param * GlobalDataBase::GetDataPara() 
{
    return globalParam;
};

void GlobalDataBase::WriteString(DataContainer *cdata, const string &cs)
{
    int nlen = static_cast<int>(cs.length());
    cdata->write(&nlen, sizeof(int));

    char *data = new char[nlen+1];

    cs.copy(data, nlen);
    data[nlen] = '\0';

    cdata->write(data, nlen+1);

    delete [] data;
}

void GlobalDataBase::ReadString(DataContainer *cdata, string &cs)
{
    int nlen = 0;
    cdata->read(&nlen, sizeof(int));

    char *data = new char[nlen+1];
    cdata->read(data, nlen+1);

    cs = data;

    delete [] data;
}

void GlobalDataBase::WriteVoid(DataContainer *cdata, void *data, int type, int size)
{
    if (type == PHSTRING)
    {
        string * mydata = reinterpret_cast<string *> (data);
        for (int i = 0; i < size; ++ i)
        {
            WriteString(cdata, mydata[i]);
        }
    }
    else if (type == PHINT)
    {
        int *mydata = reinterpret_cast<int *> (data);
        for (int i = 0; i < size; ++ i)
        {
            cdata->write(&mydata[i], sizeof(int));
        }
    }
    else if (type == PHFLOAT)
    {
        RFloat *mydata = reinterpret_cast<RFloat *> (data);
        for (int i = 0; i < size; ++ i)
        {
            cdata->write(&mydata[i], sizeof(RFloat));
        }
    }
    else if (type == PHDOUBLE)
    {
        RDouble *mydata = reinterpret_cast<RDouble *> (data);
        for (int i = 0; i < size; ++ i)
        {
            cdata->write(&mydata[i], sizeof(RDouble));
        }
    }
}

void GlobalDataBase::UpdateData(DataContainer *cdata, const string &name, int type, int size)
{
    if (type == PHSTRING)
    {
        string *data = new string [size];
        for (int i = 0; i < size; ++ i)
        {
            ReadString(cdata, data[i]);
        }
        UpdateData(name, data, type, size);
        delete [] data;
    }
    else if (type == PHINT)
    {
        int *data = new int [size];
        for (int i = 0; i < size; ++ i)
        {
            cdata->read(&data[i], sizeof(int));
        }
        UpdateData(name, data, type, size);
        delete [] data;
    }
    else if (type == PHFLOAT)
    {
        RFloat *data = new RFloat [size];
        for (int i = 0; i < size; ++ i)
        {
            cdata->read(&data[i], sizeof(RFloat));
        }
        UpdateData(name, data, type, size);
        delete [] data;
    }
    else if (type == PHDOUBLE)
    {
        RDouble *data = new RDouble [size];
        for (int i = 0; i < size; ++ i)
        {
            cdata->read(&data[i], sizeof(RDouble));
        }
        UpdateData(name, data, type, size);
        delete [] data;
    }
}

void GlobalDataBase::CompressData(DataContainer *&cdata)
{
    set < Data_SafeData > *basep = GetBaseP();
    set < Data_SafeData >::iterator it;

    uint_t ndata = basep->size();

    cdata->write(&ndata, sizeof(uint_t));

    for (it = basep->begin(); it != basep->end(); ++ it)
    {
        int type = it->GetType();
        int size = it->GetSize();

        WriteString(cdata, it->GetName());
        cdata->write(&type, sizeof(int));
        cdata->write(&size, sizeof(int));
        WriteVoid(cdata, it->GetData(), type, size);
    }
}

void GlobalDataBase::DecompressData(DataContainer *cdata)
{
    //! It seems that move cdata to the beginning is very important.
    cdata->moveto_begin();

    uint_t ndata = 0;
    cdata->read(&ndata, sizeof(uint_t));

    for (uint_t i = 0; i < ndata; ++ i)
    {
        string name;
        int type;
        int size;

        ReadString(cdata, name);
        cdata->read(&type, sizeof(int));
        cdata->read(&size, sizeof(int));
        UpdateData(cdata, name, type, size);
    }
}

string GlobalDataBase::GetStrParaFromDB(const string &name)
{
    string data;
    GetData(name, &data, PHSTRING, 1);
    return data;
}

int GlobalDataBase::GetIntParaFromDB(const string &name)
{
    int data;
    GetData(name, &data, PHINT, 1);
    return data;
}

RDouble GlobalDataBase::GetDoubleParaFromDB(const string &name)
{
    RDouble data;
    GetData(name, &data, PHDOUBLE, 1);
    return data;
}

bool GlobalDataBase::IsExist(const string &name, int type, int size)
{
    return globalParam->IsExist(name, type, size);
}

RDouble * GlobalDataBase::GetRDoubleArrayFromDB(const string &name, int numberOfElements)
{
    RDouble *rRDoubleArray = new RDouble[numberOfElements];
    GetData(name, rRDoubleArray, PHDOUBLE, numberOfElements);
    return rRDoubleArray;
}

void GlobalDataBase::GetRDoubleVectorFromDB(vector< RDouble > &parameter, const string &nameOfRFloatArray, int numberOfElements)
{
    parameter.resize(numberOfElements);
    GetData(nameOfRFloatArray, &parameter[0], PHDOUBLE, numberOfElements);
}

LIB_EXPORT int GetTaskCode()
{
    int nsimutask = GlobalDataBase::GetIntParaFromDB("nsimutask");
    return nsimutask;
}

}