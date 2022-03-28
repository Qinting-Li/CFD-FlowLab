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
//! @file      Pre_CGNSToPlot3D_Struct.h
//! @brief     Grid conversion from CGNS to Plot3D.
//! @author    Refactored by He Xianyao, Xu Gang, He Xin.

#pragma once
#include "LIB_Macro.h"
#include "Precision.h"
#include "cgnslib.h"
#include <set>

using namespace std;

namespace PHSPACE
{

class CGNS_Str_Data;

//! @brief Pre_CGNSToPlot3D_Struct class defines the method of Grid conversion from CGNS format to Plot3D.
//! Pre_CGNSToPlot3D_Struct
//! CGNS_Str_Data, BaseData
class Pre_CGNSToPlot3D_Struct
{
private:
    //! Grid file name of original grid that need to be converted.
    string gridFileName;

public:
    //! @param[in] gridFileName    Grid file name of CGNS grid.
    //! mode = 0 (output .grd, .inp files); mode(Default) = 1 (output .grd, .inp, .bc, .bcname files with the name "_0").
    LIB_EXPORT Pre_CGNSToPlot3D_Struct(const string &gridFileName, int mode = 1);

    LIB_EXPORT ~Pre_CGNSToPlot3D_Struct();

public:
    //! Run to convert CGNS grid.
    LIB_EXPORT void Run();

private:
    //! Read the grid information from CGNS file.
    void ReadCGNSGrid();

    //! Write the plot3d format grid file.
    void WriteToGrd();

    //! Write the boundary condition file.
    void WriteToInp();

    //! Write the boundary name file.
    void WriteToBcname();

    //! Write the boundary face information into file.
    void WriteToBcinformation();

    //! 
    void ResetGridScaleAndTranslate();

    //!
    void writetoCGNS(fstream &bcfile, cgsize_t *pnts, int lab, string bcName, int nCoords);

    //!
    void tram(int *a, int *b, int *v, int nCoords);

    //!
    int del(int x, int y);

private:
    CGNS_Str_Data *Str_CGNS_Data;

    string cel_file;
    string bnd_file;
    string bcname_file;
    string bcinfor_file;

    int nzones;

    set<pair<int, string> > bcNamePairSet;

    int m_mode;
};

class BaseData;
typedef char char_33[33];
typedef cgsize_t cgsize_t_60000[60000];

//! @brief Collection the grid information from CGNS.
//! CGNS_Str_Data
//! BaseData
class CGNS_Str_Data
{
public:
    CGNS_Str_Data(int nzones_in);
    ~CGNS_Str_Data(void);

private:
    int nzones;
    BaseData **base_cgns;

public:
    BaseData * GetBaseData(int izone) { return base_cgns[izone]; }
};

//! @brief Collection the grid information in each block.
//! BaseData
class BaseData
{
public:
    BaseData(void);
    ~BaseData(void);

private:
    int nCoords, idim, jdim, kdim;
    int nconnect, nBCRegions, n1to1;
    RDouble ***x, ***y, ***z;
    char_33 *zonenames, *familyName;
    BCType_t *bocoType;
    cgsize_t_60000 *pnts;
    string *bcName;
    int *lab1, *lab2;
    int ier;
    cgsize_t **srange, **donor_range;
    int **b, **vt;
    char_33 *connectname, *donorname;

public:
    void SetIdim(int data_in) { this->idim = data_in; }
    void SetJdim(int data_in) { this->jdim = data_in; }
    void SetKdim(int data_in) { this->kdim = data_in; }

    int GetIdim() { return idim; }
    int GetJdim() { return jdim; }
    int GetKdim() { return kdim; }

    void SetNCoords(int nCoords_in) { this->nCoords = nCoords_in; }
    int GetNCoords() { return nCoords; }

    void SetX(RDouble ***x_in) { this->x = x_in; }
    void SetY(RDouble ***y_in) { this->y = y_in; }
    void SetZ(RDouble ***z_in) { this->z = z_in; }
    RDouble *** GetX() { return x; }
    RDouble *** GetY() { return y; }
    RDouble *** GetZ() { return z; }

    void SetZonenames(char_33 *zonenames_in) { this->zonenames = zonenames_in; }
    char_33 * GetZonenames() { return zonenames; }

    void SetFamilyName(char_33 *familyName_in) { this->familyName = familyName_in; }
    char_33 * GetFamilyName() { return familyName; }

    void SetBocoType(BCType_t *bocoType_in) { this->bocoType = bocoType_in; }
    BCType_t * GetBocoType() { return bocoType; }

    void SetPnts(cgsize_t_60000 *pnts_in) { this->pnts = pnts_in; }
    cgsize_t_60000 * GetPnts() { return pnts; }

    void SetBcName(string *bcName_in) { this->bcName = bcName_in; }
    string * GetBcName() { return bcName; }

    void SetLab1(int *lab1_in) { this->lab1 = lab1_in; }
    int * GetLab1() { return lab1; }

    void SetNconnect(int nconnect_in) { this->nconnect = nconnect_in; }
    int GetNconnect() { return nconnect; }

    void SetNBCRegions(int nBCRegions_in) { this->nBCRegions = nBCRegions_in; }
    int GetNBCRegions() { return nBCRegions; }

    void SetN1to1(int n1to1_in) { this->n1to1 = n1to1_in; }
    int GetN1to1() { return n1to1; }

    void SetB(int *b_in, int one21_in) { this->b[one21_in] = b_in; }
    int * GetB(int one21_in) { return b[one21_in]; }

    void SetVt(int *vt_in, int one21_in) { this->vt[one21_in] = vt_in; }
    int * GetVt(int one21_in) { return vt[one21_in]; }

    void ResizeRange(int n1to1_in);
    void SetSrange(cgsize_t *srange_in, int one21_in) { this->srange[one21_in] = srange_in; }
    cgsize_t * GetSrange(int one21_in) { return srange[one21_in]; }

    void SetDonor_range(cgsize_t *donor_range_in, int one21_in) { this->donor_range[one21_in] = donor_range_in; }
    cgsize_t * GetDonor_range(int one21_in) { return donor_range[one21_in]; }

    void SetLab2(int lab2_in, int one21_in) { this->lab2[one21_in] = lab2_in; }
    int Getlab2(int one21_in) { return lab2[one21_in]; }

    void SetIer(int ier_in) { this->ier = ier_in; }
    int GetIer() { return ier; }

    void SetConnectname(char_33 *connectname_in) { this->connectname = connectname_in; }
    void SetDonorname(char_33 *donorname_in) { this->donorname = donorname_in; }
    char_33 * GetConnectname() { return connectname; }
    char_33 * GetDonorname() { return donorname; }
};

}