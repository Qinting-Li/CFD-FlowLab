#include "DataContainer.h"
using namespace std;

namespace PHSPACE
{

DataChar::DataChar()
{
    data = new vector<char>;
    pos  = 0;
}

DataChar::~DataChar()
{
    delete data;
}

void DataChar::moveto(CharVecSizeType pos)
{
    if ((0 < pos) && (pos < size()))
    {
        this->pos = pos;
    }
    else if (0 == pos)
    {
        this->pos = pos;
    }
    else
    {
        ostringstream oss;
        oss << "Out of Range: pos " << pos << endl;

        //PHMPI::FreeBlockData();
        cout << oss.str() << endl;
        exit(0);
    }
}

void DataChar::forward_position(CharVecSizeType size)
{
    pos += size;
}

CharVecSizeType DataChar::size()
{
    return data->size();
}

char * DataChar::begin()
{ 
    //! should consider that the size = 0;
    if (this->size() == 0) return 0;
    return &((*data)[0]);
}

char * DataChar::data_pointer()
{ 
    //! should consider that the size = 0;
    if (this->size() == 0) return 0;
    return &((*data)[pos]);
}

void DataChar::to_string(string &str)
{
    if (size())
    {
        str.append(begin(), size());
    }
}

void DataChar::write(void *data, CharVecSizeType size)
{
    if (size <= 0) return;

    memcpy(data_pointer(), data, size);
    forward_position(size);
}

void DataChar::read(void *data, CharVecSizeType size)
{
    if (size <= 0) return;
    memcpy(data, data_pointer(), size);
    forward_position(size);
}

void DataChar::write(void *data, CharVecSizeType size, int pos)
{
    moveto(pos);
    write(data, size);
}

void DataChar::read(void *data, CharVecSizeType size, int pos)
{
    moveto(pos);
    read(data, size);
}

void DataChar::resize(CharVecSizeType size)
{
    this->data->resize(size);
}

void DataChar::readfile(fstream &file)
{
    CharVecSizeType nlen = this->size();
    if (nlen <= 0) return;

    char *data = this->begin();

    file.read(data, nlen);
}

void DataChar::writefile(fstream &file)
{
    CharVecSizeType nlen = this->size();
    if (nlen <= 0) return;

    char *data = this->begin();

    file.write(data, nlen);
}

DataContainer::DataContainer(DataContainer *rhs)
{
    max_unit_size = 102400000;
    start  = 0;
    finish = 1;
    pos = 0;
    elem_pos = 0;
    data = new vector< DataChar * >;
    ContainerSizeType nelement = rhs->element_size();
    data->resize(nelement);
    for (ContainerSizeType i = 0; i < nelement; ++ i)
    {
        DataChar *rhsData = rhs->getiter(i);
        CharVecSizeType nlen = rhsData->size();

        DataChar *sectionData = new DataChar;
        sectionData->resize(nlen);

        sectionData->moveto_begin();
        rhsData->moveto_begin();
        memcpy(sectionData->data_pointer(), rhsData->data_pointer(), nlen);
        //sectionData->forward_position(nlen);
        (*data)[i] = sectionData;

        //this->forward_position(nlen);
    }
}

DataContainer::DataContainer()
{
    //max_unit_size = 1024000;
    max_unit_size = 102400000;
    start  = 0;
    finish = 1;
    data = new vector< DataChar * >;
    data->push_back(new DataChar());
    pos = 0;
    elem_pos = 0;
}

DataContainer::~DataContainer()
{
    for (ContainerSizeType i = 0; i < data->size(); ++ i)
    {
        delete (*data)[i];
    }

    delete data;
}

ContainerSizeType DataContainer::element_size()
{
    return data->size();
}

DataChar * DataContainer::getcurr()
{
    //! concern that the current element is self consistent with the pos.
    elem_pos = pos / max_unit_size;
    return (*data)[elem_pos];
}

DataChar * DataContainer::getiter(ContainerSizeType i)
{
    return (*data)[i];
}

void DataContainer::forward_position(ContainerSizeType size)
{
    pos += size;
}

void DataContainer::write(DataContainer *dataContainer)
{
    dataContainer->moveto_begin();
    CharVecSizeType length = dataContainer->size();

    char *data;

    data = new char[length];

    dataContainer->read(data, length);

    this->write(data, length);

    delete [] data;
}

void DataContainer::read(void *data, CharVecSizeType size)
{
    //! note. reading from DataContainer require the DataContainer is not empty.
    if (size <= 0) return;

    CharVecSizeType remain_size = remain();

    if (remain_size >= size)
    {
        //! if the sub data length is not small than 'size', it can read directly.
        getcurr()->read(data, size);
        forward_position(size);
    }
    else
    {
        //! read the data by segmentation.
        //! first read the remained data int this DataChar.
        getcurr()->read(data, remain_size);
        forward_position(remain_size);

        void *new_data = move_pointer(data, remain_size);
        CharVecSizeType new_size = size - remain_size;

        read(new_data, new_size);
    }
}

void DataContainer::write(void *data, CharVecSizeType size)
{
    if (size <= 0) return;

    secure_relative_space(size);

    CharVecSizeType remain_size = remain();

    if (remain_size >= size)
    {
        getcurr()->write(data, size);
        forward_position(size);
    }
    else
    {
        //! write the data by segmentation.
        //! write some data for the rest space.
        getcurr()->write(data, remain_size);

        forward_position(remain_size);

        void *new_data = move_pointer(data, remain_size);
        CharVecSizeType new_size = size - remain_size;

        //! repeat this function for the rest data that has node been written.
        write(new_data, new_size);
    }
}

void DataContainer::read_string(string &cs)
{
    int nlen = 0;
    read(&nlen, sizeof(int));

    char *data = new char[nlen+1];

    read(data, nlen+1);

    cs = data;

    delete [] data;
}

void DataContainer::write_string(string &cs)
{
    int nlen = static_cast<int>(cs.length());

    write(&nlen, sizeof(int));

    char *data = new char[nlen+1];

    cs.copy(data, nlen);
    data[nlen] = '\0';

    write(data, nlen+1);

    delete [] data;
}

void DataContainer::append_string(string &cs)
{
    moveto_end();
    write_string(cs);
}

void DataContainer::write(ostringstream *oss)
{
    string str = oss->str();
    int str_size = static_cast<int>(str.size());
    write(const_cast<char *>(str.c_str()), str_size * sizeof(char));
}

CharVecSizeType DataContainer::size()
{
    CharVecSizeType sum = 0;
    for (unsigned int i = 0; i < element_size(); ++ i)
    {
        sum += getiter(i)->size();
    }
    return sum;
}

void DataContainer::resize(CharVecSizeType nlen)
{
    if (nlen <= 0) return;

    //23 divided by 3 is 7, remainder 2.

    ContainerSizeType basic_elem = nlen / max_unit_size;
    ContainerSizeType remainder  = nlen % max_unit_size;

    ContainerSizeType additional_elem = 0;
    if (remainder) additional_elem = 1;

    ContainerSizeType new_size = basic_elem + additional_elem;
    element_resize(new_size);

    for (ContainerSizeType i = 0; i < element_size(); ++ i)
    {
        CharVecSizeType need_size = max_unit_size;
        if (i == basic_elem)
        {
            need_size = remainder;
        }
        getiter(i)->resize(need_size);
    }
}

void DataContainer::element_resize(ContainerSizeType new_size)
{
    if (new_size <= element_size())
    {
        erase(new_size, element_size());
        data->resize(new_size);
    }
    else
    {
        ContainerSizeType ist = element_size();
        ContainerSizeType ied = new_size;

        for (ContainerSizeType i = ist; i != ied; ++ i)
        {
            data->push_back(new DataChar());
        }
    }
}

void DataContainer::secure_relative_space(CharVecSizeType size)
{
    CharVecSizeType need_size = pos + size;

    secure_absolute_space(need_size);
}

void DataContainer::secure_absolute_space(CharVecSizeType need_size)
{
    //! if the space is enough, then it does not need to resize.
    //! it may cause some questions of these chars, if we do not resize£¬there may be some extra chars not been cleared in the memory.
    //if (need_size <= size()) return;

    resize(need_size);
}

void DataContainer::moveto_begin()
{
    pos = 0;
    for (ContainerSizeType i = 0; i < element_size(); ++ i)
    {
        getiter(i)->moveto_begin();
    }
}

void DataContainer::moveto_end()
{
    pos = size();
    for (ContainerSizeType i = 0; i < element_size(); ++ i)
    {
        getiter(i)->moveto_end();
    }
}

CharVecSizeType DataContainer::remain()
{
    CharVecSizeType remainder  = pos % max_unit_size;
    return max_unit_size - remainder;
}

void DataContainer::readfile(fstream &file)
{
    //! write the data to DataContainer from the file.
    streamsize nlen = 0;
    file.read(reinterpret_cast<char *>(&nlen), sizeof(streamsize));
    if (nlen <= 0) return;

    secure_absolute_space(static_cast< CharVecSizeType >(nlen));

    //cout << "element_size = " << element_size() << "\n";

    for (ContainerSizeType i = 0; i < element_size(); ++ i)
    {
        getiter(i)->readfile(file);
    }
}

void DataContainer::writefile(fstream &file)
{
    streamsize nlen = size();

    //! it need to be written to the file whenever the nlen is less than 0 or not.
    file.write(reinterpret_cast<char *>(&nlen), sizeof(streamsize));
    if (nlen <= 0) return;

    //cout << "element_size = " << element_size() << "\n";

    for (ContainerSizeType i = 0; i < element_size(); ++ i)
    {
        getiter(i)->writefile(file);
    }
}

void DataContainer::to_string(string &str)
{
    for (ContainerSizeType i = 0; i < element_size(); ++ i)
    {
        //cout << " i = " << i << " to_string element_size = " << element_size() << "\n";
        getiter(i)->to_string(str);
    }
}

void DataContainer::append(void *data, CharVecSizeType size)
{
    moveto_end();
    CharVecSizeType need_size = this->size() + size;
    secure_absolute_space(need_size);

    getcurr()->write(data, size);
}

void DataContainer::destroy(DataChar *iter)
{
    delete iter;
}

void DataContainer::erase(ContainerSizeType ist, ContainerSizeType ied)
{
    for (ContainerSizeType i = ist; i != ied; ++ i)
    {
        destroy(getiter(i));
    }
}

//start to add by Xu Qingxin 2013.03.20., modified 2013.05.09.
CharVecSizeType DataContainer::my_size(streamsize size)
{
    CharVecSizeType sum = 0;
    ContainerSizeType elem_beg = static_cast< ContainerSizeType >(pos / max_unit_size);
    ContainerSizeType elem_end = static_cast< ContainerSizeType >((pos + size - 1) / max_unit_size);

    for (ContainerSizeType i = elem_beg; i <= elem_end; ++ i)
    {
        sum += getiter(i)->size();
    }

    return sum;
}

void DataContainer::init_buffer(CharVecSizeType buf_beg, int elem_beg, int elem_num)
{
    pos = elem_beg * max_unit_size + buf_beg;

    for (int i = elem_beg; i < elem_beg + elem_num; ++ i)
    {
        getiter(i)->moveto_begin();
    }
}

void DataContainer::init_buffer(int elem_num, int ptr_dim, 
     ContainerSizeType *ptr_elem_beg, ContainerSizeType *ptr_elem_len, 
     ContainerSizeType *ptr_elem_num)
{
    element_resize(elem_num);

    for (int i = 0; i < ptr_dim; i ++)
    {
        CharVecSizeType remainder  = ptr_elem_len[i] % max_unit_size;
        ContainerSizeType elem_end = ptr_elem_beg[i] + ptr_elem_num[i];

        for (ContainerSizeType inum = ptr_elem_beg[i]; inum < elem_end; inum ++)
        {
            CharVecSizeType need_size;
            if (inum == elem_end - 1)
            {
                need_size = remainder;
            }
            else
            {
                need_size = max_unit_size;
            }

            getiter(inum)->resize(need_size);
        }
    }

    moveto_begin();
}

void DataContainer::static_write(CharVecSizeType begin, void *data, CharVecSizeType size)
{
    //using namespace PHMPI;

    pos += begin;

    //check:
    if (size > this->my_size(size))
    {
        ostringstream oss;

        oss << "Error: Proc ? " //<< GetCurrentProcessorID() 
            << ", The write buffer size not enough!\n"
            << "   ---- begin =" << pos << ", actual size = " << size 
            << ", buffer size = " << this->my_size(size) << " ----\n"
            << " (In static_write, from ActionKey.cpp)\n";

        string errmsg = oss.str();

        //PHMPI::FreeBlockData();
        cout << errmsg << endl;
        exit(0);
    }

    //begin                  = begin%max_unit_size;
    begin = pos % max_unit_size;
    CharVecSizeType remain_size = remain();

    if (remain_size >= size)
    {
        getcurr()->write(data, size, static_cast<int>(begin));
        forward_position(size);
    }
    else 
    {
        //! write the data by segmentation.
        //! write some data for the rest space.
        getcurr()->write(data, remain_size, static_cast<int>(begin));

        forward_position(remain_size);

        void *new_data = move_pointer(data, remain_size);
        CharVecSizeType new_size = size - remain_size;

        //! repeat this function for the rest data that has node been written.
        static_write(new_data, new_size);
    }
}

void DataContainer::static_write(void *data, CharVecSizeType size)
{
    //using namespace PHMPI;

    //check:
    if (size > this->my_size(size))
    {
        ostringstream oss;

        oss << "Error: Proc ? " //<< GetCurrentProcessorID() 
            << ", The write buffer size not enough!\n"
            << "   ---- pos =" << pos << ", actual size = " << size 
            << ", buffer size = " << this->my_size(size) << " ----\n"
            << " (In static_write, from ActionKey.cpp)\n";

        string errmsg = oss.str();

        //PHMPI::FreeBlockData();
        cout << errmsg << endl;
        exit(0);
    }

    CharVecSizeType remain_size = remain();

    if (remain_size >= size)
    {
        getcurr()->write(data, size);
        forward_position(size);
    }
    else 
    {
        //! write the data by segmentation.
        //! write some data for the rest space.
        getcurr()->write(data, remain_size);

        forward_position(remain_size);

        void *new_data = move_pointer(data, remain_size);
        CharVecSizeType new_size = size - remain_size;

        //! repeat this function for the rest data that has node been written.
        static_write(new_data, new_size);
    }
}

void DataContainer::static_read(CharVecSizeType begin, void *data, CharVecSizeType size)
{
    //using namespace PHMPI;

    pos += begin;

    //check:
    if (size > this->my_size(size))
    {
        ostringstream oss;

        oss << "Error: Proc ? " //<< GetCurrentProcessorID() 
            << ", The read buffer size not enough!\n"
            << "   ---- begin =" << pos << ", actual size = " << size 
            << ", buffer size = " << this->my_size(size) << " ----\n"
            << " (In static_read, from ActionKey.cpp)\n";

        string errmsg = oss.str();

        //PHMPI::FreeBlockData();
        cout << errmsg << endl;
        exit(0);
    }

    CharVecSizeType remain_size = remain();
    //begin                  = begin%max_unit_size;
    begin = pos % max_unit_size;

    if (remain_size >= size)
    {
        //! if the sub data length is not small than 'size', it can read directly.
        getcurr()->read(data, size, static_cast<int>(begin));
        forward_position(size);
    }
    else
    {
        //! read the data by segmentation.
        //! first read the remained data int this superchar.
        getcurr()->read(data, remain_size, static_cast<int>(begin));
        forward_position(remain_size);

        void *new_data = move_pointer(data, remain_size);
        CharVecSizeType new_size = size - remain_size;

        static_read(new_data, new_size);
    }
}

void DataContainer::static_read(void *data, CharVecSizeType size)
{
    //using namespace PHMPI;

    streamsize need_size = pos + size;

    //check:
    if (need_size > static_cast<streamsize>(this->size()))
    {
        ostringstream oss;

        oss << "Error: Proc ? " //<< GetCurrentProcessorID() 
            << ", The read buffer size not enough!\n"
            << "   ---- pos =" << pos << ", actual size = " << need_size 
            << ", buffer size = " << this->size() << " ----\n"
            << " (In static_read, from ActionKey.cpp)\n";

        string errmsg = oss.str();

        //PHMPI::FreeBlockData();
        cout << errmsg << endl;
        exit(0);
    }

    CharVecSizeType remain_size = remain();

    if (remain_size >= size)
    {
        //! if the sub data legth is not small than 'size', it can read directly.
        getcurr()->read(data, size);
        forward_position(size);
    }
    else
    {
        //! read the data by segmentation.
        //! first read the remained data int this superchar.
        getcurr()->read(data, remain_size);
        forward_position(remain_size);

        void *new_data = move_pointer(data, remain_size);
        CharVecSizeType new_size = size - remain_size;

        static_read(new_data, new_size);
    }
}

char * move_pointer(void *data, streamsize size)
{
    return reinterpret_cast<char *>(data) + size;
}
}
