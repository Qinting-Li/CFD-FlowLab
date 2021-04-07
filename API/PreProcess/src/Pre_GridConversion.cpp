#include "BC_Unstruct.h"
#include "IO_FileName.h"
#include "Pre_HDF5File.h"
#include "Glb_Dimension.h"
#include "LinkStruct.h"
#include "Post_Visual.h"
#include "Pre_GridConversion.h"
#include "HyList.h"
#include "PHIO.h"
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{
LIB_EXPORT Pre_GridConversion::Pre_GridConversion(const string & gridFileName)
{
	this->gridFileName = gridFileName;
	grids = 0;
    nBlocks = 0;
}

LIB_EXPORT Pre_GridConversion::~Pre_GridConversion()
{
	
}

LIB_EXPORT void Pre_GridConversion::Clear()
{
	for ( int izone = 0; izone < nBlocks; ++ izone )
	{
		delete grids[izone];
	}
	delete [] grids;
}

LIB_EXPORT void Pre_GridConversion::Run()
{
	ReadGrid();

	Conversion();
}

void Pre_GridConversion::WriteFaceBoundaryName(const string & out_grid_file, Grid ** grids_in)
{
    PrintToWindow("    Writing face boundary name to .bcname file ....\n");
	Grid ** grids = grids_in;
	string from_gfile = "grid.fts";
	GlobalDataBase::GetData("from_gfile",&from_gfile,PHSTRING,1);

	set< pair<int, string> > bcNamePairSet;
	//int nBlocks = 1;
	for ( int iZone = 0; iZone < nBlocks; ++ iZone )
	{
		UnstructGrid * grid = UnstructGridCast(grids[0]);
		int nBoundFace = grid->GetNBoundFace();		

		BCRecord ** bcr = grid->GetBCRecord();
		for ( cgsize_t iFace = 0; iFace < nBoundFace; iFace ++ )
		{
			int bcType    = bcr[iFace]->GetKey();
			const string & bcName = bcr[iFace]->GetBCName();
			pair<int, string> bcNamePair(bcType, bcName);
			bcNamePairSet.insert(bcNamePair);
		}

	}

    DumpFaceBoundaryName(bcNamePairSet, gridFileName);
}

void Pre_GridConversion::WriteCellToNode(const string & targetGridFileName, Grid ** grids_in)
{
	fstream file;

	Grid ** grids = grids_in;

	string cellToNodeFileName = ChangeExtensionOfFileName( targetGridFileName, "c2n" );

	PHSPACE::OpenSerialFile( file, cellToNodeFileName, ios_base::out|ios_base::binary );

	int nBlocks = 1;
	PHWrite( file, nBlocks );

	VirtualFile * virtualFile = new VirtualFile( &file );
	virtualFile->BeginWriteWork();

	for ( int iZone = 0; iZone < nBlocks; ++ iZone )
	{
		PHWrite( virtualFile, iZone );

		int nTotalCell = grids[0]->GetNTotalCell();
		PHWrite( virtualFile, nTotalCell );

		int * cell2Node = UnstructGridCast(grids[0])->GetCell2Node();
		int * node_number_of_each_cell = UnstructGridCast(grids[0])->GetNodeNumberOfEachCell();

		int nPoints;
		int count = 0;
		for ( int icell = 0; icell < nTotalCell; icell ++ )
		{
			nPoints = node_number_of_each_cell[icell];
			PHWrite( virtualFile, nPoints );

			for ( int jPoint = 0; jPoint < nPoints; ++ jPoint )
			{
				int pointIndex = cell2Node[count];
				PHWrite( virtualFile, pointIndex );
				count ++;
			}
		}
	}

	virtualFile->EndWriteWork();
	delete virtualFile;

	PHSPACE::CloseFile( file );	
}

void Pre_GridConversion::WriteFaceBC(const string &out_grid_file, Grid ** grids_in)
{
	fstream file;

	Grid ** grids = grids_in;

	string BCFileName = ChangeExtensionOfFileName( out_grid_file, "bc" );

	PHSPACE::OpenSerialFile( file, BCFileName, ios_base::out|ios_base::binary );

	//int nBlocks = 1;
	PHWrite( file, nBlocks );

	for ( int iZone = 0; iZone < nBlocks; ++ iZone )
	{
		UnstructGrid * grid = UnstructGridCast(grids[iZone]);

		grid->WriteFaceBC(file);
	}

	PHSPACE::CloseFile( file );	
}

LIB_EXPORT void DumpGrid(const string &gridFileName, Grid **grids, int nBlocks)
{
    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
    if (dumpOldGrid)
    {
        DumpOldGrid(gridFileName, grids, nBlocks);
        return;
    }

    DumpHDF5Grid(gridFileName, grids, nBlocks);
}

LIB_EXPORT void DumpOldGrid(const string &gridFileName, Grid **grids, int nBlocks)
{
    PrintToWindow("    Writing old version grid to file ....\n");
    fstream file;
    OpenSerialFile(file, gridFileName.c_str(),ios_base::out|ios_base::binary|ios_base::trunc);

    using namespace PHMPI;

    SetNumberofGlobalZones(nBlocks);

    int nZones = GetNumberofGlobalZones();

    CreateZoneProcessorID(nZones);
    CreateZoneGridID(nZones);
    CreateZoneGridType(nZones);

    int *block_proc     = GetZoneProcessorID();
    int *block_idx      = GetZoneGridID();
    int *block_type     = GetZoneGridType();
    int *block_proc_inp = GetZoneProcessorID_INP();

    for (int izone = 0; izone < nZones; ++ izone)
    {
        block_proc[izone] = izone;
        block_idx [izone] = izone;
        block_type[izone] = grids[izone]->Type();
    }

    if (block_proc_inp)
    {
        for (int izone = 0; izone < nZones; ++ izone)
        {
            SetZoneProcessorID(izone, block_proc_inp[izone]);
        }
    }

    file.write(reinterpret_cast<char*>(&nBlocks),sizeof(int));

    int m_block_proc = PHMPI::GetZoneDistributionMethod();
    if(m_block_proc == DETERMIN_BY_PARTITION && IsConvertGridToMixGrid())
    {
        int * block_proc_dump     = GetZoneProcessorIDDump();
        file.write(reinterpret_cast<char*>(block_proc_dump),nZones*sizeof(int));
    }
    else
    {
        file.write(reinterpret_cast<char*>(block_proc),nZones*sizeof(int));
    }

    file.write(reinterpret_cast<char*>(block_idx ),nZones*sizeof(int));
    file.write(reinterpret_cast<char*>(block_type),nZones*sizeof(int));

    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        grids[izone]->WriteGrid(file);
    }

    file.close();
    file.clear();
}

void GetGridsBoundingBox(Grid **grids, int nBlocks, RDouble *pmin, RDouble *pmax)
{
    pmin[0] = LARGE;
    pmin[1] = LARGE;
    pmin[2] = LARGE;

    pmax[0] = - LARGE;
    pmax[1] = - LARGE;
    pmax[2] = - LARGE;

    for (int nb = 0; nb < nBlocks; ++ nb)
    {
        RDouble *local_pmin = grids[nb]->GetMinBox();
        RDouble *local_pmax = grids[nb]->GetMaxBox();

        pmin[0] = MIN(pmin[0], local_pmin[0]);
        pmin[1] = MIN(pmin[1], local_pmin[1]);
        pmin[2] = MIN(pmin[2], local_pmin[2]);

        pmax[0] = MAX(pmax[0], local_pmax[0]);
        pmax[1] = MAX(pmax[1], local_pmax[1]);
        pmax[2] = MAX(pmax[2], local_pmax[2]);
    }
}

void GetGridsMinMaxDS(Grid **grids, int nBlocks, RDouble &mindis, RDouble &maxdis)
{
    mindis =   LARGE;
    maxdis = - LARGE;

    for (int nb = 0; nb < nBlocks; ++ nb)
    {
        RDouble dismin, dismax;
        grids[nb]->GetMinMaxDS(dismin, dismax);

        mindis = MIN(mindis, dismin);
        maxdis = MAX(maxdis, dismax);
    }
}

void MatchInterface(Grid *grid, LinkStruct *link)
{
	InterfaceInfo * iinfo = grid->GetInterfaceInfo();
	if ( iinfo == 0 ) return;

	int nIFace    = iinfo->GetNIFace();
	int zoneIndex = grid->GetZoneID();
	int npeoridic = 0;

	VVInt &zoneid  = link->GetZoneID();
	VVInt &faceid  = link->GetFaceID();
	VVInt &facemap = link->GetFaceMap();

	int * interFace2ZoneID      = iinfo->GetInterFace2ZoneID();
	int * interFace2InterFaceID = iinfo->GetInterFace2InterFaceID();

	for ( int iiface = 0; iiface < nIFace; ++ iiface )
	{
        //! image is interface mapped from i to the global,the interface is for all blocks,,not only the local block.
        int ipos = facemap[zoneIndex][ iiface ];
		//int find_flag = false;
		bool find_flag = false;     //! Corrected by clz, 2012.3.15.

        bool tmpout = false;

        if ( zoneid[ipos].size() != 2 )
        {
            cout << " more than or less than two faces coincide\n";
            cout << " zoneIndex= " << zoneIndex << "\n";
            cout << " zoneid[ipos].size() = " << zoneid[ipos].size() << "\n";
            cout << " iiface = " << iiface << " nIFace = " << nIFace << "\n";
            tmpout = true;
        }

		for ( std::size_t m = 0; m < zoneid[ipos].size(); ++ m )
		{
            int izone = zoneid[ipos][m];
            int iface = faceid[ipos][m];
		    if ( ( izone != zoneIndex ) || ( iface != iiface ) )
			{
				interFace2ZoneID[iiface]      = izone;
				interFace2InterFaceID[iiface] = iface;
				find_flag = true;
				break;
			}
		}

		if ( !find_flag )
		{
			//! It shows that there is  non-physical artificial periodic boundary condition in this block.
			//! Of course,there may be problems in the process.
            cout << "There is an error in connection between the calculation area of this block and the butt block " << 
                    "(There may also be non-physical artificial periodic boundary conditions), please check the grid!\n";
			++ npeoridic;
		}
	}
}

void ConstructGlobalInterfaceLink(Grid **grids, int nBlocks)
{
    RDouble mindis, maxdis, pmin[3], pmax[3];

    GetGridsBoundingBox(grids, nBlocks, pmin, pmax);
    GetGridsMinMaxDS(grids, nBlocks, mindis, maxdis);

    cout << "  Min && Max edge length: " << mindis << ", " << maxdis << "\n";

    typedef DataStruct_AdtTree<int, RDouble> AdtTree;

    if (GetDim() == TWO_D)
    {
        pmin[2] -= mindis;
        pmax[2] += maxdis;
    }

    LinkStruct link(pmin, pmax, mindis / 2.0, nBlocks);

    int omit_no_bound_bc = 1;
    GlobalDataBase::GetData("omit_no_bound_bc", &omit_no_bound_bc, PHINT, 1);

    if (!omit_no_bound_bc)
    {
        //! Covert NO_BOUNDARY_CONDITION to INTERFACE.
        for (int nb = 0; nb < nBlocks; ++ nb)
        {
            grids[nb]->ChangeBCType(PHENGLEI::NO_BOUNDARY_CONDITION, PHENGLEI::INTERFACE);

            //! Consider the periodic boundary as the interface.
            grids[nb]->ChangeBCType(PHENGLEI::PERIODIC, PHENGLEI::INTERFACE);
		}
	}

    //! Build the initial link information.
    for (int nb = 0; nb < nBlocks; ++ nb)
    {
        grids[nb]->BuildGridLink(&link);
    }

    //! Further process.
    for (int nb = 0; nb < nBlocks; ++ nb)
    {
        MatchInterface(grids[nb], &link);
    }
}

void DumpCharacteristicBoundary(const string &gridFileName, Grid **grids, int nBlocks)
{
    string filename = ChangeExtensionOfFileName( gridFileName, "bcmesh" );
    filename = AddSymbolToFileName(filename, "_", 0);
    vector < DataContainer * > datalist;
    datalist.resize(0);
    

    for ( int izone =  0; izone < nBlocks; ++ izone )
    {
        DataContainer *data = new DataContainer();
        data->moveto_begin();
        
        Grid *grid = grids[izone];
        ComputeCharacteristicBoundary(grid, data);

        datalist.push_back(data);
    }

    bool CharacteristicBoundary = true;

    DumpToVTK DumpToVTK(datalist, filename, CharacteristicBoundary);
    DumpToVTK.Run();

    //for ( int izone =  0; izone < nBlocks; ++ izone )
    //{      
    //    delete datalist.at(izone);
    //}
    //datalist.resize(0);

}

void ComputeCharacteristicBoundary(Grid *grid_in, DataContainer *data)
{
    int GridID = grid_in->GetGridID()->GetIndex();
    PHWrite(data, GridID);

    int dimension = grid_in->GetDim();
    PHWrite(data, dimension);

    int type = grid_in->Type();
    PHWrite(data, &type, 1);

    if (type == UNSTRUCTGRID)
    {
        ComputeCharacteristicBoundaryUns(grid_in, data);
    }
    else if (type == STRUCTGRID)
    {
        ComputeCharacteristicBoundaryStr(grid_in, data);
    }

    int TecioMission = TEC_SPACE::Writecomplete;
    PHWrite(data, &TecioMission, 1);
}

void ComputeCharacteristicBoundaryStr(Grid *grid_in, DataContainer *data)
{
    StructGrid * grid = StructGridCast(grid_in);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int s_st[3],s_ed[3];

    RDouble3D & xx = * grid->GetStructX();
    RDouble3D & yy = * grid->GetStructY();
    RDouble3D & zz = * grid->GetStructZ();

    CompositeBCRegion * composite_bcregion = grid->GetCompositeBCRegion();
    int nBCRegion = composite_bcregion->GetnBCRegion();

    int TecioMission;
    for ( int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion )
    {
        BCRegion * bcregion = composite_bcregion->GetBCRegion(ibcregion);
        int bctype = bcregion->GetBCType();

        if ( bctype == PHENGLEI::INTERFACE || bctype < 0  ) continue;

        TecioMission = TEC_SPACE::WriteBoundary;
        PHWrite(data, &TecioMission, 1);

        PHWrite(data, &bctype, 1);

        string bcName = bcregion->GetBCName();
        data->write_string(bcName);

        int * s_lr3d = bcregion->GetFaceDirectionIndex();
        int nDim = GetDim();

        for ( int m = 0; m < nDim; ++ m )
        {
            s_st[m] = bcregion->GetStartPoint( m );
            s_ed[m] = bcregion->GetEndPoint( m );

            if ( bcregion->GetFaceDirection() == m )
            {
                if ( s_st[m] > 1 || s_lr3d[m] == 1 )
                {
                    s_st[m] += 1;
                    s_ed[m] += 1;
                }
            }
            else
            {
                s_ed[m] += 1;
            }
        }

        int ist = s_st[0];
        int ied = s_ed[0];
        int jst = s_st[1];
        int jed = s_ed[1];
        int kst = s_st[2];
        int ked = s_ed[2];

        Range I(ist,ied);
        Range J(jst,jed);
        Range K(kst,ked);

        if ( nDim != THREE_D ) 
        //if ( nk == 1 ) 
        {
            K.setRange(1,1);
            kst = 1;
            ked = 1;
        }

        Int3D newindex(I,J,K,fortranArray);

        int nsurf = bcregion->GetFaceDirection() + 1;
        int lr    = bcregion->GetFaceLeftOrRightIndex();

        int count = 0;
        for ( int k = kst; k <= ked; ++ k )
        {
            for ( int j = jst; j <= jed; ++ j )
            {
                for ( int i = ist; i <= ied; ++ i )
                {
                    newindex(i,j,k) = count;
                    ++ count;
                }
            }
        }

        int nTotalNode = count;
        int nTotalCell = 1;

        //for ( int m = 0; m < nDim; ++ m )
        for ( int m = 0; m < 3; ++ m )
        {
            nTotalCell *= ( bcregion->GetEndPoint(m) - bcregion->GetStartPoint(m) + 1 );
        }

        PHWrite(data, &nTotalNode, 1);
        PHWrite(data, &nTotalCell, 1);


        RDouble * xBoundary = new RDouble [nTotalNode];
        RDouble * yBoundary = new RDouble [nTotalNode];
        RDouble * zBoundary = new RDouble [nTotalNode];
                
        int iNode = 0;
        for ( int k = kst; k <= ked; ++ k )
        {
            for ( int j = jst; j <= jed; ++ j )
            {
                for ( int i = ist; i <= ied; ++ i )
                {
                    xBoundary[iNode] = xx(i,j,k);
                    yBoundary[iNode] = yy(i,j,k);
                    zBoundary[iNode] = zz(i,j,k);
                    
                    ++ iNode;
                }
            }
        }

        PHWrite(data, xBoundary, nTotalNode);
        PHWrite(data, yBoundary, nTotalNode);
        PHWrite(data, zBoundary, nTotalNode);

        delete [] xBoundary;
        delete [] yBoundary;
        delete [] zBoundary;

        
        int il1,jl1,kl1;
        GetNsurfIndex(nsurf, il1, jl1, kl1);

        il1 = 1 - il1;
        jl1 = 1 - jl1;
        kl1 = 1 - kl1;

        if ( nk == 1 ) kl1 = 0;

        int * cell2node = new int [ nTotalCell * 4 ];
        int ncount = 0;
        for ( int k = kst; k <= ked - kl1; ++ k )
        {
            for ( int j = jst; j <= jed - jl1; ++ j )
            {
                for ( int i = ist; i <= ied - il1; ++ i )
                {
                    int p1,p2,p3,p4;

                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    if ( nsurf == 1 )
                    {
                        p1 = newindex(i,j  ,k ) + 1;
                        p2 = newindex(i,jl, k ) + 1;
                        p3 = newindex(i,jl, kl) + 1;
                        p4 = newindex(i,j  ,kl) + 1;
                    }
                    else if ( nsurf == 2 )
                    {
                        p1 = newindex(i ,j,k ) + 1;
                        p2 = newindex(il,j,k ) + 1;
                        p3 = newindex(il,j,kl) + 1;
                        p4 = newindex(i ,j,kl) + 1;
                    }
                    else
                    {
                        p1 = newindex(i ,j ,k) + 1;
                        p2 = newindex(il,j ,k) + 1;
                        p3 = newindex(il,jl,k) + 1;
                        p4 = newindex(i ,jl,k) + 1;
                    }

                    cell2node[ncount++] = p1;
                    cell2node[ncount++] = p2;
                    cell2node[ncount++] = p3;
                    cell2node[ncount++] = p4;
                }
            }
        }

        PHWrite(data, cell2node, nTotalCell * 4);
        delete [] cell2node;
    }

    if ( grid_in->GetDim() == PHSPACE::TWO_D )
    {
        /*int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();*/

        /*RDouble3D & xx = * grid->GetStructX();
        RDouble3D & yy = * grid->GetStructY();
        RDouble3D & zz = * grid->GetStructZ();*/

        int nTotalCell = grid->GetNTotalCell();
        int nTotalNode = grid->GetNTotalNode();

        TecioMission = TEC_SPACE::WriteBlock;
        PHWrite(data, &TecioMission, 1);

        PHWrite(data, &nTotalNode, 1);
        PHWrite(data, &nTotalCell, 1);

        RDouble * xBoundary = new RDouble [nTotalNode];
        RDouble * yBoundary = new RDouble [nTotalNode];
        RDouble * zBoundary = new RDouble [nTotalNode];
    
        int ist, ied, jst, jed, kst, ked;
        grid->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

        int iNode = 0;
        for ( int k = kst; k <= ked; ++ k )
        {
            for ( int j = jst; j <= jed; ++ j )
            {
                for ( int i = ist; i <= ied; ++ i )
                {
                    xBoundary[iNode] = xx(i,j,k);
                    yBoundary[iNode] = yy(i,j,k);
                    zBoundary[iNode] = zz(i,j,k);

                    iNode ++ ;
                }
            }
        }

        PHWrite(data, xBoundary, nTotalNode);
        PHWrite(data, yBoundary, nTotalNode);
        PHWrite(data, zBoundary, nTotalNode);

        delete [] xBoundary;
        delete [] yBoundary;
        delete [] zBoundary;

        int * cell2node = new int [ nTotalCell * 4 ];
        int ncount = 0;
        for ( int j = 1; j <= nj - 1; ++ j )
        {
            for ( int i = 1; i <= ni - 1; ++ i )
            {
                cell2node[ncount++] = i     + ( j - 1 ) * ni;
                cell2node[ncount++] = i + 1 + ( j - 1 ) * ni;
                cell2node[ncount++] = i + 1 + ( j     ) * ni;
                cell2node[ncount++] = i     + ( j     ) * ni;
            }
        }

        PHWrite(data, cell2node, nTotalCell * 4);
        delete [] cell2node;
    }
}

void ComputeCharacteristicBoundaryUns(Grid *grid_in, DataContainer *data)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);

    if ( grid_in->GetDim() == PHSPACE::TWO_D )
    {
        int TecioMission;

        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nTotalNode = grid_in->GetNTotalNode();
        int nTotalCell = grid_in->GetNTotalCell();
        int nBoundFace = grid_in->GetNBoundFace();
        int nTotalFace = grid_in->GetNTotalFace();

        int * node_number_of_each_face = grid->GetNodeNumberOfEachFace();
        int * face2node = grid->GetFace2Node();
        int * left_cell_of_face = grid->GetLeftCellOfFace();
        int * right_cell_of_face = grid->GetRightCellOfFace();

        RDouble * x = grid_in->GetX();
        RDouble * y = grid_in->GetY();
        RDouble * z = grid_in->GetZ();

        vector < vector < int > > face2nodelist(nBoundFace);
        vector < int > linkmap;

        set< pair<int, string> > bcnameMap;

        BCRecord **bcr = grid->GetBCRecord();
        for ( int iface = 0; iface < nBoundFace; ++ iface )
        {
            int bctype = bcr[iface]->GetKey();
            string bcname = bcr[iface]->GetBCName();

            bcnameMap.insert(pair<int,string>(bctype, bcname));
        }

        set< pair<int, string> >::iterator iter;
        for ( iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter )
        {
            if ( (*iter).first == PHENGLEI::INTERFACE || (*iter).first < 0  ) continue;

            TecioMission = TEC_SPACE::WriteBoundary;
            PHWrite(data, &TecioMission, 1);

            int BCtype = (*iter).first;
            PHWrite(data, &BCtype, 1);

            string bcName = (*iter).second;
            data->write_string(bcName);

            face2nodelist.resize(0);
            linkmap.resize(0);

            GetFace2NodeList(grid, *iter, linkmap, face2nodelist);
            
            int NumPts = static_cast<int>(linkmap.size());
            int NumElements = static_cast<int>(face2nodelist.size());

            PHWrite(data, &NumPts, 1);
            PHWrite(data, &NumElements, 1);


            RDouble * xBoundary = new RDouble [NumPts];
            RDouble * yBoundary = new RDouble [NumPts];
            RDouble * zBoundary = new RDouble [NumPts];

            for ( int iNode = 0; iNode < NumPts; ++ iNode )
            {
                int it = linkmap[iNode];

                xBoundary[iNode] = x[it];
                yBoundary[iNode] = y[it];
                zBoundary[iNode] = z[it];
            }

            PHWrite(data, xBoundary, NumPts);
            PHWrite(data, yBoundary, NumPts);
            PHWrite(data, zBoundary, NumPts);

            delete [] xBoundary;
            delete [] yBoundary;
            delete [] zBoundary;

            int * cell2node = new int [ NumElements * 4 ];
            int ncount = 0;
            for ( int iCell = 0; iCell < NumElements; ++ iCell )
            {
                uint_t np     = face2nodelist[iCell].size();
                int index0 = face2nodelist[iCell][0];
                for ( uint_t ip = 0; ip < np; ++ ip )
                {
                    cell2node[ncount] = face2nodelist[iCell][ip] + 1;
                    ncount ++;
                }

                for ( uint_t ip = np; ip < 4; ++ ip )
                {
                    cell2node[ncount] = index0 + 1;
                    ncount ++;
                }
            }

            PHWrite(data, cell2node, NumElements * 4);
            delete [] cell2node;
        }

        HyList < int > * cell2face = new HyList < int > [ nTotalCell ];
        for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
        {
            cell2face[ iCell ].SetAverNnode( 4 );
        }
        for ( int iFace = 0; iFace < nBoundFace; ++ iFace )
        {
            int le = left_cell_of_face[ iFace ];
            cell2face[ le ].insert( iFace );
        }
        for ( int iFace = nBoundFace; iFace < nTotalFace; ++ iFace )
        {
            int le = left_cell_of_face[ iFace ];
            int re = right_cell_of_face[ iFace ];
            cell2face[ le ].insert( iFace );
            cell2face[ re ].insert( iFace );
        }


        HyList < int > * cell2node = new HyList < int > [ nTotalCell ];
        for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
        {
            cell2node[ iCell ].SetAverNnode( 4 );
        }

        for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
        {
            int face0 = cell2face[ iCell ][ 0 ];
            int ip1 = face2node[ 2 * face0 ];
            int ip2 = face2node[ 2 * face0 + 1 ];
            cell2node[ iCell ].insert( ip1 );
            cell2node[ iCell ].insert( ip2 );

            int targetNode = ip2;
            int checkNode  = ip1;
            
            int nFaceOfCell = cell2face[ iCell ].size();
            bool isClockWise = false;
            while ( ! isClockWise )
            {			
                for ( int iface = 1; iface < nFaceOfCell; ++ iface )
                {
                    int face  = cell2face[ iCell ][ iface ];
                    int node1 = face2node[ 2 * face ];
                    int node2 = face2node[ 2 * face + 1 ];

                    if ( node1 == targetNode )
                    {
                        if( node2 != checkNode )
                        {
                            cell2node[ iCell ].insert( node2 );
                            targetNode = node2;	
                        }
                        else
                        {
                            isClockWise = true;
                            break;
                        }				
                    }
                    else if ( node2 == targetNode )
                    {
                        if( node1 != checkNode )
                        {
                            cell2node[ iCell ].insert( node1 );					
                            targetNode = node1;
                        }
                        else
                        {
                            isClockWise = true;
                            break;
                        }
                    }
                }
            }
        }

        TecioMission = TEC_SPACE::WriteBlock;
        PHWrite(data, &TecioMission, 1);
        
        PHWrite(data, &nTotalNode, 1);
        PHWrite(data, &nTotalCell, 1);

        PHWrite(data, x, nTotalNode);
        PHWrite(data, y, nTotalNode);
        PHWrite(data, z, nTotalNode);


        int * cell2node_vtk = new int [ nTotalCell * 4 ];
        int ncount = 0;
        for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
        {
            if ( cell2node[ iCell ].size() == 3 )
            {
                for ( int m = 0; m < 3; ++ m )
                {
                    cell2node_vtk[ncount] = cell2node[iCell].GetData(m) + 1;
                    ncount ++;
                }
                cell2node_vtk[ncount] = cell2node[iCell].GetData(2) + 1;
                ncount ++;
            }
            else if ( cell2node[ iCell ].size() == 4 )
            {
                for ( int m = 0; m < 4; ++ m )
                {
                    cell2node_vtk[ncount] = cell2node[iCell].GetData(m) + 1;
                    ncount ++;
                }			
            }
            else
            {
                cout << "Error: this function only support triangle and quad !\n";
                exit( 0 );
            }
        }

        PHWrite(data, cell2node_vtk, nTotalCell * 4);

        delete [] cell2node_vtk;
        delete [] cell2node;
        delete [] cell2face;
    }
    else
    {
        int nTotalCell = grid_in->GetNTotalCell();
        int nBoundFace = grid_in->GetNBoundFace();

        int * node_number_of_each_face = grid->GetNodeNumberOfEachFace();
        int * left_cell_of_face = grid->GetLeftCellOfFace();
        int * right_cell_of_face = grid->GetRightCellOfFace();
        int * face2node = grid->GetFace2Node();

        RDouble * x = grid_in->GetX();
        RDouble * y = grid_in->GetY();
        RDouble * z = grid_in->GetZ();

        using namespace PHENGLEI;
        vector < vector < int > > face2nodelist(nBoundFace);
        vector < int > linkmap;

        set< pair<int, string> > bcnameMap;

        BCRecord **bcr = grid->GetBCRecord();
        for ( int iface = 0; iface < nBoundFace; ++ iface )
        {
            int bctype = bcr[iface]->GetKey();
            string bcname = bcr[iface]->GetBCName();

            bcnameMap.insert(pair<int,string>(bctype, bcname));
        }


        int TecioMission;
        set< pair<int, string> >::iterator iter;
        for ( iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter )
        {
            if ( (*iter).first == PHENGLEI::INTERFACE || (*iter).first < 0  ) continue;

            TecioMission = TEC_SPACE::WriteBoundary;
            PHWrite(data, &TecioMission, 1);

            int BCtype = (*iter).first;
            PHWrite(data, &BCtype, 1);

            string bcName = (*iter).second;
            data->write_string(bcName);

            face2nodelist.resize(0);
            linkmap.resize(0);

            GetFace2NodeList(grid, *iter, linkmap, face2nodelist);
            
            int NumPts = static_cast<int>(linkmap.size());
            int NumElements = static_cast<int>(face2nodelist.size());

            PHWrite(data, &NumPts, 1);
            PHWrite(data, &NumElements, 1);
            

            RDouble * xBoundary = new RDouble [NumPts];
            RDouble * yBoundary = new RDouble [NumPts];
            RDouble * zBoundary = new RDouble [NumPts];

            for ( int iNode = 0; iNode < NumPts; ++ iNode )
            {
                int it = linkmap[iNode];

                xBoundary[iNode] = x[it];
                yBoundary[iNode] = y[it];
                zBoundary[iNode] = z[it];
            }

            PHWrite(data, xBoundary, NumPts);
            PHWrite(data, yBoundary, NumPts);
            PHWrite(data, zBoundary, NumPts);

            delete [] xBoundary;
            delete [] yBoundary;
            delete [] zBoundary;


            int * cell2node = new int [ NumElements * 4 ];
            int ncount = 0;
            for ( int iCell = 0; iCell < NumElements; ++ iCell )
            {
                uint_t np     = face2nodelist[iCell].size();
                int index0 = face2nodelist[iCell][0];
                for ( uint_t ip = 0; ip < np; ++ ip )
                {
                    cell2node[ncount] = face2nodelist[iCell][ip] + 1;
                    ncount ++;
                }

                for ( uint_t ip = np; ip < 4; ++ ip )
                {
                    cell2node[ncount] = index0 + 1;
                    ncount ++;
                }
            }

            PHWrite(data, cell2node, NumElements * 4);
            delete [] cell2node;
        }
    }
}

}
