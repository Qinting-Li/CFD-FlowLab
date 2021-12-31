#include "Pre_CGNSToPHengLEI_Unstruct.h"
#include "GridBase.h"
#include "Math_BasisFunction.h"
#include "Geo_UnstructBC.h"
#include "PHIO.h"
#include "Glb_Dimension.h"
#include "TK_Exit.h"
using namespace std;

namespace PHSPACE
{
LIB_EXPORT Pre_CGNSToPHengLEI_Unstruct::Pre_CGNSToPHengLEI_Unstruct(const string &gridFileName) : 
    Pre_GridConversion(gridFileName)
{
    factory_cgns = 0;
}

LIB_EXPORT Pre_CGNSToPHengLEI_Unstruct::~Pre_CGNSToPHengLEI_Unstruct()
{
    delete factory_cgns;
}

void Pre_CGNSToPHengLEI_Unstruct::ReadGrid()
{
    string &cgnsfile = gridFileName;
    cout << "    CGNS file name: " << cgnsfile << "\n";

    cgsize_t isize[3][1];
    int fileid, baseid;
    cgsize_t irmin, irmax, istart, iend;
    int nbndry, iparent_flag;
    int cellDim, physDim;
    cgsize_t nTotalNode, nTotalCell;
    int nCoords, nSections, nBCRegions, nFamilies, FamBC, nGeo;
    cgsize_t nBCElem, normListFlag;
    int normalIndex, nDataSet;
    char zonename[33], sectionname[33];
    char basename[33], coordname[33], bocoName[33];
    cgsize_t elementdatasize;
    ElementType_t  etype;
    ZoneType_t     zoneType;
    DataType_t     dataType;
    BCType_t       bocoType;
    PointSetType_t ptsetType;
    DataType_t     normDataType;
    CGNSRawCoor    *raw_coor;

    //! Open the CGNS for reading and check if the file was found.
    //if ( cg_open(cgnsfile.c_str(), MODE_READ, &fileid) != CG_OK )
    if (cg_open(cgnsfile.c_str(), CG_MODE_READ, &fileid) != CG_OK)
    {
        TK_Exit::ExceptionExit(cg_get_error());
    }

    //! Determine the of bases in the grid.
    cg_nbases(fileid, &baseid);

    //! Determine the number of FamilyBC in the grid.
    cg_nfamilies(fileid, baseid, &nFamilies);

    typedef char char_33[33];
    char_33 *FamilyName;
    char_33 FamBCName;

    BCType_t *FamilyBCType = new BCType_t [nFamilies];
    int BC=1;

    if (nFamilies > 0)
    {
        FamilyName = new char_33 [nFamilies];

        for (int Fam = 1; Fam <= nFamilies; ++ Fam)
        {
            cg_family_read(fileid, baseid, Fam, FamilyName[Fam-1], &FamBC, &nGeo);
            cg_fambc_read(fileid, baseid, Fam, BC, FamBCName, &FamilyBCType[Fam-1]);
        }
    }

    //! Read the number of zones in the grid.
    cg_nzones(fileid, baseid, &nBlocks);

    //! Check the cell and physical dimensions of the bases.
    cg_base_read(fileid, baseid, basename, &cellDim, &physDim);

    if (cellDim != PHSPACE::GetDim())
    {
        TK_Exit::ExceptionExit("Error: Dimension in CGNS file is not according with that in parameter file!");
    }

    factory_cgns = new CGNSFactory(nBlocks);
    for (int izone = 1; izone <= nBlocks; ++ izone)
    {
        //! Check the zone type.
        cg_zone_type(fileid, baseid, izone, &zoneType);

        if (zoneType == Structured)
        {
            TK_Exit::ExceptionExit("Error: the input CGNS mesh must be UNSTRUCTURED type!");
        }

        //! Determine the number of vertices and volume elements in this zone.
        cg_zone_read(fileid, baseid, izone, zonename, isize[0]);

        nTotalNode = isize[0][0];
        nTotalCell = isize[1][0];

        CGNSBase *base_cgns = factory_cgns->GetCGNSBase(izone-1);

        base_cgns->SetNTotalNode(nTotalNode);
        base_cgns->SetNTotalCell(nTotalCell);

        cout << " nTotalNode = " << nTotalNode << " nTotalCell = " << nTotalCell << "\n";

        //! Determine the number and names of the coordinates.
        cg_ncoords(fileid, baseid, izone, &nCoords);

        //! Lower range index.
        irmin = 1;
        //! Upper range index of vertices.
        irmax = nTotalNode;

        //! Temporarily assume that grid is unstructured, it will be considered later.
        RDouble *x = new RDouble[nTotalNode]();
        RDouble *y = new RDouble[nTotalNode]();
        RDouble *z = new RDouble[nTotalNode]();

        RawGrid *rawgrid = factory_cgns->GetRawGrid(izone-1);

        rawgrid->SetNTotalNode(nTotalNode);
        rawgrid->SetX(x);
        rawgrid->SetY(y);
        rawgrid->SetZ(z);

        raw_coor = new CGNSRawCoor();

        for (int icoor = 0; icoor < nCoords; ++ icoor)
        {
            cg_coord_info(fileid, baseid, izone, icoor+1, &dataType, coordname);
            raw_coor->AllocateData(icoor, nTotalNode, dataType);

            //! Read the x-, y-, z-coordinates.
            cg_coord_read(fileid, baseid, izone, coordname, dataType, &irmin, &irmax, raw_coor->GetCoor(icoor));
        }

        if (nCoords == 2)
        {
            for (int iNode = 0; iNode < nTotalNode; ++ iNode)
            {
                z[iNode] = 0.0;
            }

            raw_coor->AllocateData(2, nTotalNode, dataType);
        }

        raw_coor->SetAllData(x, y, z);

        if (nCoords == 2)
        {
            for (int iNodess = 0; iNodess < nTotalNode; ++ iNodess)
            {
                z[iNodess] = 0.0;
            }
        }

        delete raw_coor;

        //! Determine the number of sections for this zone. Note that
        //! surface elements can be stored in a volume zone, but they
        //! are NOT taken into account in the number obtained from 
        //! cg_zone_read.
        cg_nsections(fileid, baseid, izone, &nSections);

        //! Loop over the number of sections and read the element
        //! connectivities. As CGNS starts the numbering at 1 the
        //! for-loop starts at 1 as well.
        base_cgns->CreateElement(nSections);

        cgsize_t *base_cgns_istart = base_cgns->GetIndexOfStart();
        cgsize_t *base_cgns_iend = base_cgns->GetIndexOfEnd();
        int *base_cgns_elem_type = base_cgns->GetElementType();
        cgsize_t **base_cgns_conn_list = base_cgns->GetElementConnectionList();
        cgsize_t  *base_cgns_elementdatasize = base_cgns->GetElementDataSize();
        cgsize_t **base_cgns_mixed_conn_list = base_cgns->GetMixedConnectionIndex();

        cgsize_t *base_elem_elem_pt = factory_cgns->GetBaseElement()->GetElementPoint();

        cout << " nSections = " << nSections << "\n";
        for (int iSection = 1; iSection <= nSections; ++ iSection)
        {
            cg_section_read(fileid, baseid, izone, iSection, sectionname, &etype,
                            &istart, &iend, &nbndry, &iparent_flag);

            cg_ElementDataSize(fileid, baseid, izone, iSection, &elementdatasize);

            cout << "\nReading section " << iSection << " ...\n";
            cout << "   section name = " << sectionname << "\n";
            cout << "   section type = " << ElementTypeName[etype] << "\n";
            cout << "   istart,iend =  " << istart << " " << iend << "\n";

            base_cgns_istart   [iSection - 1] = istart;
            base_cgns_iend     [iSection - 1] = iend;
            base_cgns_elem_type[iSection - 1] = etype;
            base_cgns_elementdatasize[iSection - 1] = elementdatasize;

            cgsize_t nelem = iend - istart + 1;
            if (etype != MIXED)
            {
                int ptnum = base_elem_elem_pt[etype];
                elementdatasize = ptnum * nelem;
            }
            base_cgns_conn_list[iSection - 1] = new cgsize_t [elementdatasize];
            
            //! Read the connectivity. Again, the node numbering of the 
            //! connectivities start at 1. If internally a starting index 
            //! of 0 is used (typical for C-codes) 1 must be subtracted 
            //! from the connectivities read.
            cg_elements_read(fileid, baseid, izone, iSection, base_cgns_conn_list[iSection - 1], NULL);

            if (etype == MIXED)
            {
                base_cgns_mixed_conn_list[iSection - 1] = new cgsize_t [nelem];

                int iNode = 0;
                int index = 0;
                for ( ; iNode < elementdatasize; )
                {
                    int ptnum = base_elem_elem_pt[base_cgns_conn_list[iSection - 1][iNode]];
                    base_cgns_mixed_conn_list[iSection - 1][index] = iNode;

                    iNode = iNode + ptnum + 1;
                    ++ index;
                }
            }
        }

        // Determine the number of boundary conditions for this zone.
        cg_nbocos(fileid, baseid, izone, &nBCRegions);

        cout << "\nNumber of BC Regions: " << nBCRegions << endl;
        base_cgns->CreateBC(nBCRegions);

        int *base_cgns_nBCElem            = base_cgns->GetNumberOfBCElements();
        int *base_cgns_bc_elem_set_type   = base_cgns->GetBoundaryElementType();
        int *base_cgns_bc_grid_location   = base_cgns->GetBoundaryGridLocation();
        int *base_cgns_bc_type            = base_cgns->GetBCType();
        cgsize_t **base_cgns_bc_conn_list = base_cgns->GetBoundaryElementConnectionList();

        GridLocation_t igr;

        //! Loop over the number of boundary conditions.
        string *boundaryName = new string [nBCRegions]();
        base_cgns->SetBCName(boundaryName);

        for (int iBCRegion = 1; iBCRegion <= nBCRegions; ++ iBCRegion)
        {
            cout << "  Reading bcRegion " << iBCRegion << " ... " << endl;
            cg_goto(fileid, baseid, "Zone_t", izone, "ZoneBC_t", 1, "BC_t", iBCRegion, "end");
            cg_gridlocation_read(&igr);

            cout << "  Grid Location Name = " << GridLocationName[igr] << "\n";
            if (igr == FaceCenter)
            {
                cout << "GridLocation = FaceCenter means BC data refers to elements, not nodes\n";
            }
            else if (igr == Vertex)
            {
                cout << "GridLocation = Vertex means BC data refers to nodes, not elements\n";
                cout << " ============================================================ " << endl;
                cout << "    Very Important Warning: " << endl;
                cout << "        if the cgns grid file is exported from Gridgen15.18,  " << endl; 
                cout << "        make sure the V16 Format button is toggle off, because" << endl;
                cout << "        the Face-based BCs are exported as Vertex-based type  " << endl;
                cout << "        wrongly by Gridgen15.18 !!!" << endl;
                cout << " =========================================================== "  << endl;
                //igr = FaceCenter;
            }
            //! Read the info for this boundary condition.
            cg_boco_info(fileid, baseid, izone, iBCRegion, bocoName, &bocoType, &ptsetType, &nBCElem, &normalIndex,
                         &normListFlag, &normDataType, &nDataSet);

            if (bocoType == FamilySpecified || nFamilies > 0)
            {
                for (int iFam = 1; iFam <= nFamilies; ++ iFam)
                {
                    if (strcmp(bocoName, FamilyName[iFam-1]) == 0)
                    {
                        bocoType = FamilyBCType[iFam-1];
                        break;
                    }
                }
            }

            boundaryName[iBCRegion - 1] = bocoName;

            cout << "  Boundary Name = " << bocoName << "\n";
            cout << "  BCtype  = " << BCTypeName[bocoType] << "\n";
            cout << "  nBCElem = " << nBCElem << "\n";

            base_cgns_nBCElem         [iBCRegion-1] = static_cast<int>(nBCElem);
            base_cgns_bc_conn_list    [iBCRegion-1] = new cgsize_t [nBCElem];
            base_cgns_bc_type         [iBCRegion-1] = bocoType;
            base_cgns_bc_elem_set_type[iBCRegion-1] = ptsetType;
            base_cgns_bc_grid_location[iBCRegion-1] = igr;

            cout << "  PointSetType_t = " << PointSetTypeName[ptsetType] << "\n";

            //! Read the element ID¡¯s.
            cg_boco_read(fileid, baseid, izone, iBCRegion, base_cgns_bc_conn_list[iBCRegion-1], NULL);

            if (nBCElem == 2)
            {
                cout << "  Start & End point of Region " << ": " 
                     << base_cgns_bc_conn_list[iBCRegion-1][0] << " " << base_cgns_bc_conn_list[iBCRegion-1][1] << "\n";
            }
            else
            {
                cout << "  The first 20 points of Region " << ": ";
                for (int j = 0; j < MIN(20, (int)nBCElem); ++ j)
                    //for ( int j = 0; j < nBCElem; ++ j )
                {
                    cout << base_cgns_bc_conn_list[iBCRegion-1][j] << " ";
                }
                cout << "\n";
            }
            //! And much more to make it fit into the 
            //! internal data structures.

            cout << endl;
        }
    }
     delete [] FamilyBCType;    FamilyBCType = NULL;
    //! Close CGNS file.
    cg_close(fileid);
}

void Pre_CGNSToPHengLEI_Unstruct::Conversion()
{
    factory_cgns->ConvertGrid2Fantasy();

    grids = new Grid * [nBlocks];

    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        GridID *index = new GridID(izone);
        Grid *grid = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
        grids[izone] = grid;
        CGNS2UnsGrid(factory_cgns, izone, UnstructGridCast(grid));
        grid->ComputeMinMaxBox();
    }
}

LIB_EXPORT void Pre_CGNSToPHengLEI_Unstruct::WriteAdditionalInformation(const string &targetGridFileName, Grid **grids_in)
{
    WriteFaceBoundaryName(targetGridFileName, grids_in);

    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
    if (dumpOldGrid)
    {
        WriteCellToNode(targetGridFileName, grids_in);

        WriteFaceBC(targetGridFileName, grids_in);
    }
}

}
