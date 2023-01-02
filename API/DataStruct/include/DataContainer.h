//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      DataContainer.h
//! @brief     It defines kernel data structure of PHengLEI, named 'DataContainer'.
//!            All types of data are standardized into DataChar. "Three-In-One" data is designed:
//!            1 - MPI communication; 2 - File I/O; 3 - Communication for hybrid parallel.
//! @author    Bell, Xu Qingxin, He Xin.

#pragma once
#include <fstream>
#include <string>
#include "TypeDefine.h"
#include "AMRDef.h"
using namespace std;

namespace PHSPACE
{
class DataChar;
typedef vector< char >::size_type CharVecSizeType;
typedef vector< DataChar * >::size_type ContainerSizeType;
typedef string::size_type StringSizeType;
//! the purpose of designing DataContainer is to make up the deficiencies of the original DataChar.
//! the interface should look like to be the same as the DataChar.
//! the design concept is as this at least.
class DataContainer
{
public:
    typedef int iterator;
public:
    DataContainer(DataContainer *rhs);    // Bell add, but i am not sure it's perfect right.
    DataContainer();
    ~DataContainer();
private:
    vector< DataChar * > *data;
    iterator start;
    iterator finish;
    CharVecSizeType pos, elem_pos, max_unit_size;

public:

    //! get data from vector *data.
    DataChar * getcurr();

    //! return vector *data's value at the site of i.
    DataChar * getiter(ContainerSizeType i);

    //! return start.
    iterator begin() { return start; }

    //! return finish.
    iterator end() { return finish; };

    //! invoking delete function to delete *iter.
    void destroy(DataChar *iter);

    //!invoking destroy function to delete *iter data,start with ist and end with ied.
    //!@param[in ]: ist  begin position.
    //!@param[in ]: ied  end position
    void erase(ContainerSizeType ist, ContainerSizeType ied);

    //! return the number of vector *data 
    ContainerSizeType element_size();

    //!read DataContainer vector<char> * data,and receive it by void *data.
    //!@param[in ]: *data , target object.
    //!@param[in ]: size
    void read(void *data, CharVecSizeType size);

    //!write *data into DataContainer vector<char> *data
    //!@param[in ]: *data ,source object.
    //!@param[in ]: size
    void write(void *data, CharVecSizeType size);

    //!initialize dataContainer with char 
    void write(DataContainer *dataContainer);

    //!read data from file.
    void readfile(fstream &file);

    //!write data into file.
    void writefile(fstream &file);

    //!read DataContainer vector<char> *data, and put it in cs.
    void read_string(string &cs);

    //!write cs into DataContainer vector<char> *data.
    void write_string(string &cs);

    //!put oss data into vector<char> *data.
    void write(ostringstream *oss);

    //!return size of sum of all data object.
    CharVecSizeType size();

    //!vector<char> *data object resize itself.
    void resize(CharVecSizeType nlen);

    //!char transform to string
    void to_string(string &str);

    //!link another data
    void append(void *data, CharVecSizeType size);

    //!link sting object
    void append_string(string &cs);

    //!before invoking resize function, calculate it size.
    void secure_relative_space(CharVecSizeType size);

    //! invoking resize function.
    void secure_absolute_space(CharVecSizeType need_size);

    //! invoking moveto function.
    void moveto_begin();

    //! invoking moveto function.
    void moveto_end();

    //! calculate: pos += size;
    void forward_position(ContainerSizeType size);

private:

    //!use new DataChar object fill up ,until the size equal new_size.
    void element_resize(ContainerSizeType new_size);

    //!return max_unit_size - remainder;
    CharVecSizeType remain();

public:

    //!invoking DataChar::read function
    //! @param[in ]:   *data :  target object.\n
    //! @param[in ]:
    void static_read(void *data, CharVecSizeType size);

    //!invoking DataChar::write function
    //! @param[in ]: *data :source object.\n
    //! @param[in ]:
    void static_write(void *data, CharVecSizeType size);


    //!invoking DataChar::read function with two parameter
    //! @param[in ]:   *data :  target object.\n
    //! @param[in ]:
    //! @param[in ]:
    void static_read(CharVecSizeType begin, void *data, CharVecSizeType size);

    //!invoking DataChar::write function with two parameter
    //! @param[in ]:  *data  :source object.\n
    //! @param[in ]:
    //! @param[in ]:
    void static_write(CharVecSizeType begin, void *data, CharVecSizeType size);
public:

    //!return size of sum of all data object
    CharVecSizeType my_size(streamsize size);

    //!
    streamsize get_remainder(streamsize position) 
        const { return position % max_unit_size; }

    //!
    int get_element_index(streamsize position) 
        const { return (int)(position / max_unit_size); }

    //!
    streamsize get_element_pos(streamsize elem) 
        const { return elem * max_unit_size; }
public:

    //!from elem_beg to elem_num.use  moveto(0) function initial vector<char> * data .
    //! @param[in ]: buf_beg   : \n
    //! @param[in ]: elem_beg  : the number of invoking moveto_begin() function.\n
    //! @param[in ]: elem_num  : the number of invoking moveto_begin() function.
    void init_buffer(CharVecSizeType buf_beg, int elem_beg, int elem_num);

    //!
    //! @param[in ]: ptr_elem_len  : use to resize vector<char> *data.\n
    //! @param[in ]: ptr_elem_beg  : use to store which object to be resize.\n
    //! @param[in ]: ptr_elem_num  : use to store which object to be resize.\n
    //! @param[in ]: elem_num
    //! @param[in ]: ptr_dim
    void init_buffer(int elem_num, int ptr_dim, ContainerSizeType *ptr_elem_beg, 
         ContainerSizeType *ptr_elem_len, ContainerSizeType *ptr_elem_num);
};

class DataChar
{
public:
    DataChar();
    ~DataChar();
public:

    //!
    CharVecSizeType size();

    //! invoking function memcpy
    //! @param[in ]:   *data : target object.\n
    //! @param[in ]:        
    void read(void *data, CharVecSizeType size);

    //! invoking function read with two parameter
    //! @param[in ]:   *data :  target object.\n
    //! @param[in ]:   pos:     moveto(pos);\n
    //! @param[in ]:
    void read(void *data, CharVecSizeType size, int pos);

    //! invoking function memcpy
    //! @param[in ]:   *data :source object.\n
    //! @param[in ]:            
    void write(void *data, CharVecSizeType size);

    //! invoking function write with two parameter
    //! @param[in ]:   *data:    source object.\n
    //! @param[in ]:   pos:      moveto(pos);\n
    //! @param[in ]:
    void write(void *data, CharVecSizeType size, int pos);

    //!file.read(data, nlen);
    void readfile(fstream &file);

    //!file.write(data, nlen);
    void writefile(fstream &file);

    //!str.append(begin(),size());
    void to_string(string &str);

    //!return &((*data)[0]);
    char * begin();

    //!
    void moveto_begin() { moveto(0); }

    //!
    void moveto_end() { pos = size(); }

    //!invoking data->resize(size);
    void resize(CharVecSizeType size);

    //!return &((*data)[pos]);
    char * data_pointer();

    //! calculate pos += size;
    void forward_position(CharVecSizeType size);

private:

    //!
    void moveto(CharVecSizeType pos);

private:
    CharVecSizeType pos;
    vector<char> *data;

public:

    //!
    char * data_pointer(int beg) { return &((*data)[beg]); }
};

//!return reinterpret_cast<char *>(data) + size;
char * move_pointer(void *data, streamsize size);
}
