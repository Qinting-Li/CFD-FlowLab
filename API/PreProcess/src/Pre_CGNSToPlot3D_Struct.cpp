#include "Pre_CGNSToPlot3D_Struct.h"
#include "Pre_GridConversion.h"
#include "PHIO.h"
#include "Pointer.h"
#include "TK_Parse.h"
#include "GridType.h"
#include "Pre_GridBase.h"
#include "TK_Exit.h"
#include "Math_BasisFunction.h"
#include "Constants.h"

namespace PHSPACE
{

Pre_CGNSToPlot3D_Struct::Pre_CGNSToPlot3D_Struct(const string &gridFileName, int mode)
{
    this->gridFileName = gridFileName;
    Str_CGNS_Data = 0;

    m_mode = mode;
    nzones = 0;
}

Pre_CGNSToPlot3D_Struct::~Pre_CGNSToPlot3D_Struct()
{
    delete Str_CGNS_Data;
}

void Pre_CGNSToPlot3D_Struct::Run()
{
    ReadCGNSGrid();

    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
    if (abs(gridScaleFactor - 1.0) > SMALL)
    {
        ResetGridScaleAndTranslate();
    }
    if (m_mode == 0) // Output .grd, .inp file with the same name.
    {
        WriteToGrd();

        WriteToInp();
    }
    else if (m_mode == 1) // Output .grd, .inp, .bcname, .bc file with the name include "_0". 
    {
        WriteToGrd();

        WriteToInp();

        int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
        if (dumpOldGrid)
        {
            WriteToBcinformation();
        }
    }
}

void Pre_CGNSToPlot3D_Struct::ReadCGNSGrid()
{
    int fileid = 0, baseid = 0, nbases = 0;
    int CellDim = 0, PhyDim = 0;
    int nfamililies = 0, nFambc = 0, nGeo = 0;
    int normalIndex[3] = {0}, nDataSet = 0;
    char_33 basename, coordname, bocoName;
    char_33 *familynames, famBCName;
    cgsize_t isize[3][3], nPoints, normListFlag;
    BCType_t *fambctype;
    PointSetType_t ptsetType;
    ZoneType_t zone_type;
    DataType_t dataType;
    DataType_t normDataType;
    int ier = 0;

    string &cgnsfile = gridFileName;
    cout << "    CGNS file name: " << cgnsfile << "\n";

    string fname = " ", fext = " ";
    GetNameExt(cgnsfile, fname, fext, ".");
    if (m_mode == 0)
    {
        cel_file = fname + ".grd";
        bnd_file = fname + ".inp";
    }
    else if (m_mode == 1) // The file name include "_0".
    {
        cel_file     = fname + "_0.grd";
        bnd_file     = fname + "_0.inp";
        bcname_file  = fname + ".bcname";
        bcinfor_file = fname + "_0.bc";

        string wallTemperatureFileName = fname + ".temperature";
        GlobalDataBase::UpdateData("wallTemperatureFileName", &wallTemperatureFileName, PHSTRING, 1);
    }

    //! update from_gfile to covert from grd to fts.
    GlobalDataBase::UpdateData("from_gfile", &cel_file, PHSTRING, 1);
    GlobalDataBase::UpdateData("bnd_file",   &bnd_file, PHSTRING, 1);

    //! Open the CGNS for reading and check if the file was found.
    if (cg_open(cgnsfile.c_str(), CG_MODE_READ, &fileid) != CG_OK)
    {
        TK_Exit::ExceptionExit(cg_get_error());
    }

    //! Determine the number of bases in the grid
    cg_nbases(fileid, &nbases);

    if (nbases > 1) 
    {
        TK_Exit::ExceptionExit("Cannot deal this situation! nbase > 1\n");
    }

    baseid = 1;
    //! Check the cell and physical dimensions of the bases.
    cg_base_read(fileid, baseid, basename, &CellDim, &PhyDim);

    //! Determine the number of FamilyBC in the grid
    cg_nfamilies(fileid, baseid, &nfamililies);
    fambctype   = new BCType_t[nfamililies];
    familynames = new char_33 [nfamililies];

    if (nfamililies > 0)
    {
        for (int fam = 1; fam <= nfamililies; ++ fam)
        {
            fambctype[fam-1] = BCTypeNull;
            cg_family_read(fileid, baseid, fam, familynames[fam-1], &nFambc, &nGeo);
            cg_fambc_read(fileid, baseid, fam, 1, famBCName, &fambctype[fam-1]);
        }
    }

    //! Read the number of zones in the grid.
    cg_nzones(fileid, baseid, &nzones);
    cout << "	Block number = " << nzones << endl;
    Str_CGNS_Data = new CGNS_Str_Data(nzones);

    for (int iZone = 1; iZone <= nzones; ++ iZone)
    {
        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone-1);

        //! Determine the number and names of the coordinates.
        base_cgns->SetNCoords(CellDim);

        char_33 *zonenames = new char_33 [1];
        base_cgns->SetZonenames(zonenames);

        cg_zone_type(fileid, baseid, iZone, &zone_type);
        if (zone_type == Unstructured)
        {
            TK_Exit::ExceptionExit("Error: the input CGNS mesh must be STRUCTURED type!");
        }

        cg_zone_read(fileid, baseid, iZone, zonenames[0], isize[0]);

        char_33 zonename;
        cg_zone_read(fileid, baseid, iZone, zonename, isize[0]);

        //! Import i¡¢j¡¢k direction information.
        cgsize_t idim = isize[0][0];
        cgsize_t jdim = isize[0][1];
        cgsize_t kdim = 1;

        if (CellDim ==3)
        {
            kdim = isize[0][2];
        }

        base_cgns->SetIdim(static_cast<int>(idim));
        base_cgns->SetJdim(static_cast<int>(jdim));
        base_cgns->SetKdim(static_cast<int>(kdim));

        cgsize_t irmin[3] = {0}, irmax[3] = {0};

        //! lower range index
        irmin[0] = 1;
        irmin[1] = 1;
        irmin[2] = 1;

        //! upper range index of vertices
        irmax[0] = idim;
        irmax[1] = jdim;
        irmax[2] = kdim;

        RDouble ***x = NewPointer3< RDouble >(static_cast<int>(kdim), static_cast<int>(jdim), static_cast<int>(idim));
        RDouble ***y = NewPointer3< RDouble >(static_cast<int>(kdim), static_cast<int>(jdim), static_cast<int>(idim));
        RDouble ***z = NewPointer3< RDouble >(static_cast<int>(kdim), static_cast<int>(jdim), static_cast<int>(idim));
        base_cgns->SetX(x);
        base_cgns->SetY(y);
        base_cgns->SetZ(z);

        int icoor = 0;
        cg_coord_info(fileid, baseid, iZone, icoor+1,  &dataType, coordname);
        cg_coord_read(fileid, baseid, iZone, coordname, dataType, irmin, irmax, **x);
        cg_coord_info(fileid, baseid, iZone, icoor+2,  &dataType, coordname);
        cg_coord_read(fileid, baseid, iZone, coordname, dataType, irmin, irmax, **y);
        if (CellDim == 3)
        {
            cg_coord_info(fileid, baseid, iZone, icoor+3,  &dataType, coordname);
            cg_coord_read(fileid, baseid, iZone, coordname, dataType, irmin, irmax, **z);
        }

        //! Read the Connectivity information
        int nBCRegions, n1to1;
        cg_goto  (fileid, baseid, "Zone_t", iZone, "ZoneBC_t", 1, "end");
        cg_nbocos(fileid, baseid, iZone, &nBCRegions);
        base_cgns->SetNBCRegions(nBCRegions);

        cg_n1to1(fileid, baseid, iZone, &n1to1);
        char_33 *connectname = new char_33 [n1to1]();
        char_33 *donorname   = new char_33 [n1to1]();
        base_cgns->SetConnectname(connectname);
        base_cgns->SetDonorname(donorname);

        int mark = 0;
        int make_edge = 0;

        //! Determine the number of interface repeats.
        for (int one21 = 1; one21 <= n1to1; ++ one21)
        {
            cgsize_t srange[6], donor_range[6];
            int transform[3];
            cg_1to1_read(fileid, baseid, iZone, one21, connectname[one21-1], donorname[one21-1], srange, donor_range, transform);
            if (strcmp(zonename, donorname[one21-1]) == 0)
            {
                mark ++;
            }
        }

        int nInterfaceBCRegions = n1to1;
        if (CellDim == 3 && n1to1 > 2)
        {
            //! CGNS Grid with more than 2 interfaces.
            nInterfaceBCRegions = n1to1 - mark/2;
        }
        base_cgns->SetN1to1(nInterfaceBCRegions);

        int      *lab1       = new int     [nBCRegions]();
        string   *bcName     = new string  [nBCRegions]();
        char_33  *familyName = new char_33 [nBCRegions]();
        BCType_t *bocoType   = new BCType_t[nBCRegions]();
        cgsize_t_60000 *pnts = new cgsize_t_60000[nBCRegions]();

        base_cgns->SetLab1(lab1);
        base_cgns->SetBocoType(bocoType);
        base_cgns->SetFamilyName(familyName);
        base_cgns->SetPnts(pnts);
        base_cgns->SetBcName(bcName);

        //! Read boundary condition information
        for (int boco = 0; boco < nBCRegions; ++ boco)
        {
            ier = cg_goto(fileid, baseid, "Zone_t", iZone, "ZoneBC_t", 1, "BC_t", boco+1, "end");
            cg_boco_info(fileid, baseid, iZone, boco+1, bocoName, &bocoType[boco], &ptsetType, &nPoints, &normalIndex[0], 
                         &normListFlag, &normDataType, &nDataSet);
            cg_boco_read(fileid, baseid, iZone, boco+1, pnts[boco], NULL);

            int make_edge_judge = 0;
            for (int nrange = 0; nrange < 3; ++ nrange)
            {
                if ((pnts[boco][nrange] - pnts[boco][nrange+3]) == 0)
                {
                    make_edge_judge ++;
                }
            }

            if (make_edge_judge >= 2)
            {
                make_edge = make_edge + 1;
            }

            cg_famname_read(familyName[boco]);

            if (make_edge_judge >= 2)
            {
                continue;
            }

            //! Set BCtypeName;
            bcName[boco] = BCTypeName[bocoType[boco]];
            if (nfamililies > 0)
            {
                for (int fam = 1; fam <= nfamililies; ++ fam)
                {
                    if (strcmp(familyName[boco], familynames[fam - 1]) == 0)
                    {
                        bocoType[boco] = fambctype[fam-1];
                        bcName[boco] = familyName[boco];
                        break;
                    }
                }
            }

            BaseBCType *bc_convert = new BaseBCType();
            lab1[boco] = bc_convert->GetBCType(bocoType[boco]);

            if (strcmp(BCTypeName[bocoType[boco]], "BCDegenerateLine") == 0)
            {
                if (CellDim == 3)
                {
                    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
                    int iMin = pnts[boco][0] - 1;
                    int jMin = pnts[boco][1] - 1;
                    int kMin = pnts[boco][2] - 1;
                    int iMax = pnts[boco][3] - 1;
                    int jMax = pnts[boco][4] - 1;
                    int kMax = pnts[boco][5] - 1;
                    RDouble dx = abs(x[kMin][jMin][iMin] - x[kMin][jMin][iMax]) * gridScaleFactor;
                    RDouble dy = abs(y[kMin][jMin][iMin] - y[kMin][jMin][iMax]) * gridScaleFactor;
                    RDouble dz = abs(z[kMin][jMin][iMin] - z[kMin][jMin][iMax]) * gridScaleFactor;
                    if ((dx+dy+dz) < SMALL && pnts[boco][0] != pnts[boco][3])
                    {
                        lab1[boco] = 71;
                    }
                    
                    dx = abs(x[kMin][jMin][iMin] - x[kMin][jMax][iMin]) * gridScaleFactor;
                    dy = abs(y[kMin][jMin][iMin] - y[kMin][jMax][iMin]) * gridScaleFactor;
                    dz = abs(z[kMin][jMin][iMin] - z[kMin][jMax][iMin]) * gridScaleFactor;
                    if ((dx+dy+dz) < SMALL && pnts[boco][1] != pnts[boco][4])
                    {
                        lab1[boco] = 72;
                    }
                    
                    dx = abs(x[kMin][jMin][iMin] - x[kMax][jMin][iMin]) * gridScaleFactor;
                    dy = abs(y[kMin][jMin][iMin] - y[kMax][jMin][iMin]) * gridScaleFactor;
                    dz = abs(z[kMin][jMin][iMin] - z[kMax][jMin][iMin]) * gridScaleFactor;
                    if ((dx+dy+dz) < SMALL && pnts[boco][2] != pnts[boco][5]) 
                    {
                        lab1[boco] = 73;
                    }
                }
                else
                {
                    lab1[boco] = 71;
                }
            }

            pair<int, string> bcNamePair(lab1[boco], bcName[boco]);
            bcNamePairSet.insert(bcNamePair);
        }

        //! Determine the number of interface connectivity entries.
        int nconnect;
        if (CellDim == 3)
        {
            nconnect = nBCRegions + nInterfaceBCRegions - make_edge;
        }
        else
        {
            nconnect = nBCRegions + nInterfaceBCRegions;
        }
        base_cgns->SetNconnect(nconnect);

        //! Zone Connectivity
        ier = cg_goto(fileid, baseid, "Zone_t", iZone, "ZoneGridConnectivity_t", 1, "GridConnectivity1to1_t", 1, "end");
        base_cgns->SetIer(ier);

        int lab2 = 0;
        if (CellDim == 3)
        {
            base_cgns->ResizeRange(nInterfaceBCRegions - make_edge);
        }
        else
        {
            base_cgns->ResizeRange(nInterfaceBCRegions);
        }

        if (ier == CG_OK)
        {
            int count = 0;

            for (int one21 = 1; one21 <= n1to1; ++ one21)
            {
                int transform[3] = {0};

                int *b = new int [3]();
                cgsize_t *srange, *donor_range;
                srange = new cgsize_t [6]();
                donor_range = new cgsize_t [6]();

                cg_1to1_read(fileid, baseid, iZone, one21, connectname[one21-1], donorname[one21-1], srange, donor_range, transform);

                base_cgns->SetSrange(srange, count);
                base_cgns->SetDonor_range(donor_range, count);

                for (int zoneid = 1; zoneid <= nzones; ++ zoneid)
                {
                    cg_zone_read(fileid, baseid, zoneid, zonename, isize[0]);
                    if (strcmp(zonename, donorname[one21-1]) == 0)
                    {
                        lab2 = zoneid;
                        break;
                    }
                }

                if (lab2 == iZone)
                {
                    string connectname1 = connectname[one21-1];
                    if (strcmp(connectname1.substr(0,14).c_str(), "1to1InterfaceB") == 0 && CellDim == 3)
                    {
                        continue;
                    }
                }

                uint_t idcmp1;
                uint_t idcmp2;
                if (CellDim == 3)
                {
                    if (srange[0] == srange[3])
                    {
                        idcmp1 = abs(srange[1] - srange[4]);
                        idcmp2 = abs(srange[2] - srange[5]);
                        if (idcmp1 < idcmp2) 
                        {
                            b[0] =  1;
                            b[1] = -1;
                            b[2] =  1;
                        }
                        else
                        {
                            b[0] =  1;
                            b[1] =  1;
                            b[2] = -1;
                        }
                    }
                    else if (srange[1] == srange[4])
                    {
                        idcmp1 = abs(srange[2] - srange[5]);
                        idcmp2 = abs(srange[0] - srange[3]);
                        if (idcmp1 < idcmp2) 
                        {
                            b[0] =  1;
                            b[1] =  1;
                            b[2] = -1;
                        }
                        else
                        {
                            b[0] = -1;
                            b[1] =  1;
                            b[2] =  1;
                        }
                    }
                    else if (srange[2] == srange[5])
                    {
                        idcmp1 = abs(srange[0] - srange[3]);
                        idcmp2 = abs(srange[1] - srange[4]);
                        if (idcmp1 < idcmp2)
                        {
                            b[0] = -1;
                            b[1] =  1;
                            b[2] =  1;
                        }
                        else
                        {
                            b[0] =  1;
                            b[1] = -1;
                            b[2] =  1;
                        }
                    }
                }
                else
                {
                    if (srange[0] == srange[2])
                    {
                        b[0] =  1;
                        b[1] = -1;
                    }
                    else if (srange[1] == srange[3])
                    {
                        b[0] =  -1;
                        b[1] =  1;
                    }
                }

                int *vt = new int[3]();
                base_cgns->SetVt(vt, count);
                base_cgns->SetB(b, count);

                tram(transform, b, vt, CellDim);

                base_cgns->SetLab2(lab2, count);
                ++ count;

            }
        }
    }
    delete [] familynames;    familynames = NULL;
    delete [] fambctype;    fambctype = NULL;
    cg_close(fileid);
}

void Pre_CGNSToPlot3D_Struct::ResetGridScaleAndTranslate()
{
    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
    RDouble *gridTranslationVector = new RDouble[3]();
    if (GlobalDataBase::IsExist("gridTranslationVector", PHDOUBLE, 3))
    {
        GlobalDataBase::GetData("gridTranslationVector", gridTranslationVector, PHDOUBLE, 3);
    }

    for (int iZone = 1; iZone <= nzones; ++ iZone)
    {
        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone-1);

        int nCoords =base_cgns->GetNCoords();

        int idim = base_cgns->GetIdim();
        int jdim = base_cgns->GetJdim();
        int kdim = base_cgns->GetKdim();

        RDouble ***x = base_cgns->GetX();
        RDouble ***y = base_cgns->GetY();
        RDouble ***z = base_cgns->GetZ();

        for (int k = 0; k < kdim; ++ k)
        {
            for (int j = 0; j < jdim; ++ j)
            {
                for (int i = 0; i < idim; ++ i)
                {
                    x[k][j][i] *= gridScaleFactor;
                    x[k][j][i] += gridTranslationVector[0];
                }
            }
        }

        for (int k = 0; k < kdim; ++ k)
        {
            for (int j = 0; j < jdim; ++ j)
            {
                for (int i = 0; i < idim; ++ i)
                {
                    y[k][j][i] *= gridScaleFactor;
                    y[k][j][i] += gridTranslationVector[1];
                }
            }
        }

        for (int k = 0; k < kdim; ++ k)
        {
            for (int j = 0; j < jdim; ++ j)
            {
                for (int i = 0; i < idim; ++ i)
                {
                    if (nCoords == 2) z[k][j][i] = 0.0;
                    z[k][j][i] *= gridScaleFactor;
                    if (nCoords == 3) z[k][j][i] += gridTranslationVector[2];
                }
            }
        }
    }
}

void Pre_CGNSToPlot3D_Struct::WriteToGrd()
{
    //! Creating Grid file.
    fstream Gridgenfile;
    Gridgenfile.open(cel_file.c_str(), ios_base::out|ios_base::binary|ios_base::trunc);
    Gridgenfile.write(reinterpret_cast<char *>(&nzones), sizeof(int));

    //! Import the maximum value of i, j, k.
    for (int iZone = 1; iZone <= nzones; ++ iZone)
    {
        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone-1);

        int idim = base_cgns->GetIdim();
        int jdim = base_cgns->GetJdim();
        int kdim = base_cgns->GetKdim();

        Gridgenfile.write(reinterpret_cast<char *>(&idim), sizeof(int));
        Gridgenfile.write(reinterpret_cast<char *>(&jdim), sizeof(int));
        Gridgenfile.write(reinterpret_cast<char *>(&kdim), sizeof(int));
    }

    //! Import the coordinates of node in the CGNS Grid.
    for (int iZone = 1; iZone <= nzones; ++ iZone)
    {
        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone-1);

        int nCoords =base_cgns->GetNCoords();

        int idim = base_cgns->GetIdim();
        int jdim = base_cgns->GetJdim();
        int kdim = base_cgns->GetKdim();

        RDouble ***x = base_cgns->GetX();
        RDouble ***y = base_cgns->GetY();
        RDouble ***z = base_cgns->GetZ();

        for (int k = 0; k < kdim; ++ k)
        {
            for (int j = 0; j < jdim; ++ j)
            {
                for (int i = 0; i < idim; ++ i)
                {
                    Gridgenfile.write(reinterpret_cast<char *>(&x[k][j][i]), sizeof(RDouble));
                }
            }
        }

        for (int k = 0; k < kdim; ++ k)
        {
            for (int j = 0; j < jdim; ++ j)
            {
                for (int i = 0; i < idim; ++ i)
                {
                    Gridgenfile.write(reinterpret_cast<char *>(&y[k][j][i]), sizeof(RDouble));
                }
            }
        }

        for (int k = 0; k < kdim; ++ k)
        {
            for (int j = 0; j < jdim; ++ j)
            {
                for (int i = 0; i < idim; ++ i)
                {
                    if (nCoords == 2) z[k][j][i] = 0.0;
                    Gridgenfile.write(reinterpret_cast<char *>(&z[k][j][i]), sizeof(RDouble));
                }
            }
        }
    }

    Gridgenfile.close();
}

void Pre_CGNSToPlot3D_Struct::WriteToInp()
{
    //! Create Inp file.
    fstream GridgenBCfile;
    GridgenBCfile.open(bnd_file.c_str(), ios_base::out);
    GridgenBCfile << 1 << endl;
    GridgenBCfile << nzones << endl;

    for (int iZone = 1; iZone <= nzones; ++ iZone)
    {
        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone-1);

        int nCoords = base_cgns->GetNCoords();

        int idim = base_cgns->GetIdim();
        int jdim = base_cgns->GetJdim();
        int kdim = base_cgns->GetKdim();
        char_33 *zonename = base_cgns->GetZonenames();

        //! Set the minimum data length.
        int wordwidth = 8;

        GridgenBCfile << setw(wordwidth) << idim;
        GridgenBCfile << setw(wordwidth) << jdim;
        if (nCoords == 3)
        {
            GridgenBCfile << setw(wordwidth) << kdim;
        }
        GridgenBCfile << endl;
        GridgenBCfile << zonename[0] << endl;

        int nconnect = base_cgns->GetNconnect();
        GridgenBCfile <<  nconnect << endl;

        int    nBCRegions    = base_cgns->GetNBCRegions();
        int    *lab1         = base_cgns->GetLab1();
        string *bcName       = base_cgns->GetBcName();
        cgsize_t_60000 *pnts = base_cgns->GetPnts();

        for (int boco = 1; boco <= nBCRegions; ++ boco)
        {
            writetoCGNS(GridgenBCfile, pnts[boco - 1], lab1[boco - 1], bcName[boco - 1], nCoords);
        }

        int ier = base_cgns->GetIer();
        if (ier == CG_OK)
        {
            int n1to1 = base_cgns->GetN1to1();
            for (int one21 = 1; one21 <= n1to1; ++ one21)
            {
                int *b = base_cgns->GetB(one21-1);
                int *vt = base_cgns->GetVt(one21-1);
                cgsize_t *sranges = base_cgns->GetSrange(one21-1);
                cgsize_t *donor_ranges = base_cgns->GetDonor_range(one21-1);

                //! Import FarField Boundary Condition and Wall Boundary Condition.
                for (int iCoords = 0; iCoords < nCoords; ++ iCoords)
                {
                    GridgenBCfile << setw(8) << sranges[iCoords] * b[iCoords];
                    GridgenBCfile << setw(8) << sranges[iCoords + nCoords] * b[iCoords];
                }

                //! Import interface information.
                GridgenBCfile << setw(8) << -1;
                GridgenBCfile << "   " << "Connection" << endl;

                for (int iCoords = 0; iCoords < nCoords; ++ iCoords)
                {
                    GridgenBCfile << setw(8) << donor_ranges[iCoords] * vt[iCoords];
                    GridgenBCfile << setw(8) << donor_ranges[iCoords + nCoords] * vt[iCoords];
                }

                int lab2 = base_cgns->Getlab2(one21-1);

                //! Import the information of symmetrical interface.
                GridgenBCfile << setw(8) << lab2 << endl;
            }
        }
        else
        {
            cout << "Warning: " << "Zone " << (iZone-1) << " has no interface!" << endl;
        }
    }

    GridgenBCfile.close();
}

void Pre_CGNSToPlot3D_Struct::WriteToBcname()
{
    cout << "After Reading CGNS grid,\n";
    cout << "Number of Original Boundary Condition: " << bcNamePairSet.size() << "\n";
    int iBC = 0;
    for (set< pair<int, string> >::iterator iter = bcNamePairSet.begin(); iter != bcNamePairSet.end(); ++ iter)
    {
        string bcName = (*iter).second;
        if (bcName == "Connection") continue;

        cout << "  The " << iBC + 1 << " CGNS BC Name is: " << bcName << endl;
    }

    DumpBoundaryInfo(bcNamePairSet);
}

void Pre_CGNSToPlot3D_Struct::WriteToBcinformation()
{
    fstream file;
    PHSPACE::OpenFile(file, bcinfor_file, ios_base::out|ios_base::binary);

    PHWrite(file, nzones);

    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        VirtualFile *virtualFile = new VirtualFile(&file);
        virtualFile->BeginWriteWork();

        PHWrite(virtualFile, iZone);

        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone);

        int nBCRegions = base_cgns->GetNBCRegions();
        int nconnect = base_cgns->GetNconnect();
        PHWrite(virtualFile, nconnect);

        string *bcNameList = new string [nconnect];
        string *bcName = base_cgns->GetBcName();
        uint_long totalSize = 0;
        for (int boco = 0; boco < nBCRegions; ++ boco)
        {
            bcNameList[boco] = bcName[boco];

            totalSize += static_cast< uint_long > (bcNameList[boco].size());
            totalSize += 1;
        }

        int ier = base_cgns->GetIer();
        if (ier == CG_OK)
        {
            int n1to1 = base_cgns->GetN1to1();
            for (int one21 = 0; one21 < n1to1; ++ one21)
            {
                bcNameList[nBCRegions+one21] = "Connection";

                totalSize += static_cast< uint_long > (bcNameList[nBCRegions+one21].size());
                totalSize += 1;
            }
        }
        else
        {
            cout << "Warning: " << "Zone " << iZone << " has no interface!" << endl;
        }

        char *bcNameChar = new char[totalSize];
        unsigned int count = 0;
        for (int boco = 0; boco < nconnect; ++ boco)
        {
            string &bcName = bcNameList[boco];
            streamsize nameSize = bcName.size();
            for (unsigned int iChar = 0; iChar < nameSize; ++ iChar)
            {
                bcNameChar[count ++] = bcName[iChar];
            }
            bcNameChar[count ++] = '\0';
        }

        PHWrite(virtualFile, totalSize);
        PHWrite(virtualFile, bcNameChar, static_cast<int>(totalSize));

        virtualFile->EndWriteWork();
        delete virtualFile;
        delete [] bcNameList;
        delete [] bcNameChar;
    }

    PHSPACE::CloseFile(file);
}

void Pre_CGNSToPlot3D_Struct::writetoCGNS(fstream &bcfile, cgsize_t *pnts, int lab, string bcName, int nCoords)
{
    if (nCoords == 2)
    {
        bcfile << setw(8) << pnts[0];
        bcfile << setw(8) << pnts[2];
        bcfile << setw(8) << pnts[1];
        bcfile << setw(8) << pnts[3];
        bcfile << setw(8) << lab;
    }
    else
    {
        bcfile << setw(8) << pnts[0];
        bcfile << setw(8) << pnts[3];
        bcfile << setw(8) << pnts[1];
        bcfile << setw(8) << pnts[4];
        bcfile << setw(8) << pnts[2];
        bcfile << setw(8) << pnts[5];
        if (bcName.substr(0,7) == "BndCond")
        {
            bcfile << setw(8) << atoi(bcName.substr(8, 2).c_str()) - 17;
        }
        else
        {
            bcfile << setw(8) << lab;
        }
    }

    bcfile << "   " << bcName;
    bcfile << endl;
}

void Pre_CGNSToPlot3D_Struct::tram(int *a, int *b, int *v, int nCoords)
{
    int t[3][3];
    int i, j;
    for (i = 0; i < nCoords; ++ i)
    {
        v[i] = 0;
        for (j = 0; j < nCoords; ++ j)
        {
            t[i][j] = del(a[j], i+1);
            v[i] = v[i] + t[i][j] * b[j];
        }
    }

    return;
}

int Pre_CGNSToPlot3D_Struct::del(int x, int y)
{
    if (ABS(x) == ABS(y))
    {
        return 1;
    }
    return 0;
}

CGNS_Str_Data::CGNS_Str_Data(int nzones_in)
{
    this->nzones = nzones_in;
    base_cgns = new BaseData *[nzones_in];
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        base_cgns[iZone] = new BaseData();
    }
}

CGNS_Str_Data::~CGNS_Str_Data(void)
{
    delete [] base_cgns;
}

BaseData::BaseData(void)
{
    nCoords     = 0;
    idim        = 0;
    jdim        = 0;
    kdim        = 0;
    nconnect    = 0;
    nBCRegions  = 0;
    n1to1       = 0;
    pnts        = 0;
    zonenames   = 0;
    bcName      = 0;
    bocoType    = 0;
    familyName  = 0;
    lab2        = 0;
    srange      = 0;
    donor_range = 0;
    b           = 0;
    vt          = 0;
}

BaseData::~BaseData(void)
{
    DelPointer3(x);
    DelPointer3(y);
    DelPointer3(z);
    DelPointer2(b);
    DelPointer2(vt);
    DelPointer2(srange);
    DelPointer2(donor_range);

    delete [] bcName;
    delete [] zonenames;
    delete [] familyName;
    delete [] bocoType;
    delete [] pnts;
    delete [] lab1;
    delete [] lab2;
    delete [] connectname;
    delete [] donorname;
}

void BaseData::ResizeRange(int n1to1_in)
{
    b           = new int *[n1to1_in];
    vt          = new int *[n1to1_in];
    lab2        = new int  [n1to1_in];
    srange      = new cgsize_t *[n1to1_in];
    donor_range = new cgsize_t *[n1to1_in];
}

}