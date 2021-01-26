#include <list>
#include "MeshAgglomeration.h"
#include "GridType.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "Geo_FaceMetrics_Unstruct.h"
#include "Geo_CellMetrics_Unstruct.h"
#include "Geo_DynamicGridMetrics_Unstruct.h"
#include "Geo_LSQWeight_Unstruct.h"
#include "BC_Unstruct.h"
#include "Geo_Face.h"
#include "GRHash.h"
#include "Geo_FaceTopo_Unstruct.h"

namespace PHSPACE{
LIB_EXPORT Mesh_Agglomeration::Mesh_Agglomeration(UnstructGrid * fineGridIn, int fineGridLevelIn)
{
	this->fineGrid = fineGridIn;
	this->fineGridLevel = fineGridLevelIn;

	this->coarseGrid = 0;
}

LIB_EXPORT Mesh_Agglomeration::~Mesh_Agglomeration()
{
}

LIB_EXPORT void Mesh_Agglomeration::CoarsenGridOnce()
{
	int nTotalFace = fineGrid->GetNTotalFace();
    int nBoundFace = fineGrid->GetNBoundFace();
    int nTotalCell = fineGrid->GetNTotalCell();
    int n = nBoundFace + nTotalCell;
    int *left_cell_of_face = fineGrid->GetLeftCellOfFace();
    int *right_cell_of_face = fineGrid->GetRightCellOfFace();

    RDouble * vol = fineGrid->GetCellVolume();
    RDouble * ns  = fineGrid->GetFaceArea();

	RDouble * xcc = fineGrid->GetCellCenterX();
	RDouble * ycc = fineGrid->GetCellCenterY();
	RDouble * zcc = fineGrid->GetCellCenterZ();

	RDouble * virtualCellCenterX = new RDouble [nTotalFace];
	RDouble * virtualCellCenterY = new RDouble [nTotalFace];
	RDouble * virtualCellCenterZ = new RDouble [nTotalFace];
	PHSPACE::SetField(virtualCellCenterX, 0.0, nTotalFace);
	PHSPACE::SetField(virtualCellCenterY, 0.0, nTotalFace);
	PHSPACE::SetField(virtualCellCenterZ, 0.0, nTotalFace);


    RDouble *tarea = new RDouble[nTotalCell];   // total area of each cell
    for ( int icell = 0; icell < nTotalCell; ++ icell )
	{
		tarea[icell] = 0.0;
	}

    int count = 0;
    int le, re;
    for ( int iface = 0; iface < nBoundFace; ++ iface )
    {
        le = left_cell_of_face[ iface ];
        tarea[le] += ns[ iface ];
    }

    for ( int iface = nBoundFace; iface < nTotalFace; ++ iface )
    {
        le = left_cell_of_face [ iface ];
        re = right_cell_of_face[ iface ];
        tarea [le] += ns[ iface ];
        tarea [re] += ns[ iface ];
    }

	int * fmark = new int[nTotalFace];

	SetField(fmark, -1, nTotalFace);
	int cnTFace = 0;

	//! Prepare the inner-faces that have same left&&right cell, on coarse grid.
	//! fmark: mapping face on fine grid to face on coarse grid.
	//! fmark = -1: the two neighbor cells would not been agglomerated.
	using namespace GRIDER_SPACE;
	int nInteriorFace = nTotalFace - nBoundFace;
	FaceHash faceHash;
	vector<Geo_Face *> faceList;
	for ( int iface = nBoundFace; iface < nTotalFace; ++ iface )
	{
		le = left_cell_of_face [ iface ];
		re = right_cell_of_face[ iface ];

		Geo_Face * face = new Geo_Face(le, re, iface, nInteriorFace);

		bool isFaceExist;

		// Insert the node into the tetNode, which is a hash table.
		int faceIndexinHash = faceHash.insert(face, isFaceExist) - 1;

		// Push the node into node list of GRAdvanceFront.
		if(isFaceExist)
		{
			// the coarse face have been built.
			delete face;
			face = faceHash[faceIndexinHash];
			fmark[iface] = face->GetID();
		}else
		{
			// build new coarse face.
			face->SetID( cnTFace );
			fmark[iface] = cnTFace;
			cnTFace ++;

			faceList.push_back(face);
		}
	}
	faceHash.clear();
	vector<Geo_Face *>::iterator faceIter;
	for (faceIter = faceList.begin(); faceIter != faceList.end(); ++ faceIter)
	{
		delete (* faceIter);
	}
	faceList.clear();


	// barea: total area between two cells, may be more than one face between two cells on coarse grid.
    RDouble *barea = new RDouble[ cnTFace ];
    for ( int i = 0; i < cnTFace; ++ i )
	{
		barea[i] = 0.0;
	}

    for ( int iface = nBoundFace; iface < nTotalFace; ++ iface )
    {
        barea[ fmark[iface] ] += ns[iface];
    }

    RDouble * ratio = new RDouble[nTotalFace];
    int fno;
    for ( int iface = nBoundFace; iface < nTotalFace; ++ iface )
    {
        le = left_cell_of_face [ iface ];
        re = right_cell_of_face[ iface ];
        fno = fmark[ iface ];

        ratio[iface] = pow( ( tarea[le] + tarea[re] - 2.0 * barea[fno] ), 1.5 ) / ( vol[le] + vol[re] );

		virtualCellCenterX[iface] = (vol[le] * xcc[le] + vol[re] * xcc[re]) / ( vol[le] + vol[re] );
		virtualCellCenterY[iface] = (vol[le] * ycc[le] + vol[re] * ycc[re]) / ( vol[le] + vol[re] );
		virtualCellCenterZ[iface] = (vol[le] * zcc[le] + vol[re] * zcc[re]) / ( vol[le] + vol[re] );
    }


    delete [] barea;
    delete [] tarea;

	RDouble limitAngle = -90.0;
// 	if(fineGridLevel >= 1)
// 	{
// 		limitAngle = -50.0 * fineGridLevel;
// 	}
	AspectRatioModifyBySkewness(ratio, virtualCellCenterX, virtualCellCenterY, virtualCellCenterZ, limitAngle);

    typedef pair<RDouble,int> TRatioIndex;
    typedef set<TRatioIndex,less<TRatioIndex> > Tset;

    Tset myset;

    for ( int iface = nBoundFace; iface < nTotalFace; ++ iface )
    {
        myset.insert( TRatioIndex(ratio[iface], iface) );
    }

    Tset::iterator curr = myset.begin();

    typedef list<int> Tlist;
    Tlist mylist;

    bool * unMerged = new bool [nTotalCell];
    for ( int icell = 0; icell < nTotalCell; ++ icell )
	{
		unMerged[icell] = true;
	}

    for ( curr = myset.begin(); curr != myset.end(); ++ curr )
    {
        fno = (*curr).second;
        le = left_cell_of_face[ fno ];
        re = right_cell_of_face[ fno ];
		if ( unMerged[le] && unMerged[re] && ratio[fno] >= 0.0)
        {
            mylist.push_back(fno);
            unMerged[le] = false;
            unMerged[re] = false;
        }
    }

	delete [] virtualCellCenterX;
	delete [] virtualCellCenterY;
	delete [] virtualCellCenterZ;

    // now merging
    int *cell2coarsegridcell = new int[ nTotalCell ];
    fineGrid->SetCell2CoarseGridCell(cell2coarsegridcell);
    for ( int icell = 0; icell < nTotalCell; ++ icell )
	{
		cell2coarsegridcell[icell] = - 1;
	}

    int cnTCell = 0;
    Tlist::iterator lcurr;
    for ( lcurr = mylist.begin(); lcurr != mylist.end(); ++ lcurr )
    {
        fno =  * lcurr;
        le = left_cell_of_face[ fno ];
        re = right_cell_of_face[ fno ];
        cell2coarsegridcell[le] = cnTCell;
        cell2coarsegridcell[re] = cnTCell ++;
    }

    // some cells may not have been merged, merge them now
	for ( int icell = 0; icell < nTotalCell; ++ icell )
	{
		if(cell2coarsegridcell[icell] == -1)
		{
			cell2coarsegridcell[icell] = cnTCell ++;
		}		
	}

    delete [] ratio;
    delete [] unMerged;
	delete [] fmark;

	coarseGrid = new UnstructGrid();
	coarseGrid->InitGrid(new GridID(*fineGrid->GetGridID()), fineGridLevel + 1, fineGrid->GetDim(), UNSTRUCTGRID);

    coarseGrid->SetNTotalCell( cnTCell );
    int nTotalNode = fineGrid->GetNTotalNode();
    coarseGrid->SetNTotalNode( nTotalNode );
    RDouble *x, *y, *z;
    RDouble * xx = fineGrid->GetX();
    RDouble * yy = fineGrid->GetY();
    RDouble * zz = fineGrid->GetZ();

    x = new RDouble[ nTotalNode ];
    y = new RDouble[ nTotalNode ];
    z = new RDouble[ nTotalNode ];

    for ( int inode = 0; inode < nTotalNode; ++ inode )
    {
        x[inode] = xx[inode];
        y[inode] = yy[inode];
        z[inode] = zz[inode];
    }
    coarseGrid->SetX(x);
    coarseGrid->SetY(y);
    coarseGrid->SetZ(z);
}

//global agglomeration method used to coarsen the grid to generate cell2coarsegridcell only
LIB_EXPORT void Mesh_Agglomeration::CoarsenGridOnceWithoutSkewnessLimit()
{
	int nTotalFace = fineGrid->GetNTotalFace();
	int nBoundFace = fineGrid->GetNBoundFace();
	int nTotalCell = fineGrid->GetNTotalCell();
	int n = nBoundFace + nTotalCell;
	int *left_cell_of_face = fineGrid->GetLeftCellOfFace();
	int *right_cell_of_face = fineGrid->GetRightCellOfFace();

	RDouble * vol = fineGrid->GetCellVolume();
	RDouble * ns  = fineGrid->GetFaceArea();

	RDouble * xcc = fineGrid->GetCellCenterX();
	RDouble * ycc = fineGrid->GetCellCenterY();
	RDouble * zcc = fineGrid->GetCellCenterZ();

	int count = 0;
	int le, re;

	RDouble *tarea = new RDouble[nTotalCell];   // total area of each cell
	for ( int icell = 0; icell < nTotalCell; ++ icell )
	{
		tarea[icell] = 0.0;
	}

	for ( int iface = 0; iface < nBoundFace; ++ iface )
	{
		le = left_cell_of_face[ iface ];
		tarea[le] += ns[ iface ];
	}

	for ( int iface = nBoundFace; iface < nTotalFace; ++ iface )
	{
		le = left_cell_of_face [ iface ];
		re = right_cell_of_face[ iface ];
		tarea [le] += ns[ iface ];
		tarea [re] += ns[ iface ];
	}

	int * fmark = new int[nTotalFace];
	SetField(fmark, -1, nTotalFace);
	int cnTFace = 0;

	//! Prepare the inner-faces that have same left&&right cell, on coarse grid.
	//! fmark: mapping face on fine grid to face on coarse grid.
	//! fmark = -1: the two neighbor cells would not been agglomerated.
	using namespace GRIDER_SPACE;
	int nInteriorFace = nTotalFace - nBoundFace;
	FaceHash faceHash;
	vector<Geo_Face *> faceList;
	for ( int iface = nBoundFace; iface < nTotalFace; ++ iface )
	{
		le = left_cell_of_face [ iface ];
		re = right_cell_of_face[ iface ];

		Geo_Face * face = new Geo_Face(le, re, iface, nInteriorFace);

		bool isFaceExist;

		// Insert the node into the tetNode, which is a hash table.
		int faceIndexinHash = faceHash.insert(face, isFaceExist) - 1;

		// Push the node into node list of GRAdvanceFront.
		if(isFaceExist)
		{
			// the coarse face have been built.
			delete face;
			face = faceHash[faceIndexinHash];
			fmark[iface] = face->GetID();
		}else
		{
			// build new coarse face.
			face->SetID( cnTFace );
			fmark[iface] = cnTFace;
			cnTFace ++;
		}
	}
	faceHash.clear();
	vector<Geo_Face *>::iterator faceIter;
	for (faceIter = faceList.begin(); faceIter != faceList.end(); ++ faceIter)
	{
		delete (* faceIter);
	}
	faceList.clear();

	// barea: total area between two cells, may be more than one face between two cells on coarse grid.
	RDouble *barea = new RDouble[ cnTFace ];
	for ( int i = 0; i < cnTFace; ++ i )
	{
		barea[i] = 0.0;
	}
	for ( int iface = nBoundFace; iface < nTotalFace; ++ iface )
	{
		barea[ fmark[iface] ] += ns[iface];
	}

	RDouble * ratio = new RDouble[nTotalFace];
	int fno;
	for ( int iface = nBoundFace; iface < nTotalFace; ++ iface )
	{
		le = left_cell_of_face [ iface ];
		re = right_cell_of_face[ iface ];
		fno = fmark[ iface ];
		ratio[iface] = pow( ( tarea[le] + tarea[re] - 2.0 * barea[fno] ), 1.5 ) / ( vol[le] + vol[re] );
	}
	delete [] barea;
	delete [] tarea;


	typedef pair<RDouble,int> TRatioIndex;
	typedef set<TRatioIndex,less<TRatioIndex> > Tset;
	Tset myset;
	for ( int iface = nBoundFace; iface < nTotalFace; ++ iface )
	{
		myset.insert( TRatioIndex(ratio[iface], iface) );
	}

	Tset::iterator curr = myset.begin();
	typedef list<int> Tlist;
	Tlist mylist;
	bool * unMerged = new bool [nTotalCell];
	for ( int icell = 0; icell < nTotalCell; ++ icell )
	{
		unMerged[icell] = true;
	}
	for ( curr = myset.begin(); curr != myset.end(); ++ curr )
	{
		fno = (*curr).second;
		le = left_cell_of_face[ fno ];
		re = right_cell_of_face[ fno ];
		if ( unMerged[le] && unMerged[re])
		{
			mylist.push_back(fno);
			unMerged[le] = false;
			unMerged[re] = false;
		}
	}

	// now merging
	int *cell2coarsegridcell = new int[ nTotalCell ];
	fineGrid->SetCell2CoarseGridCell(cell2coarsegridcell);
	for ( int icell = 0; icell < nTotalCell; ++ icell )
	{
		cell2coarsegridcell[icell] = - 1;
	}

	int cnTCell = 0;
	Tlist::iterator lcurr;
	for ( lcurr = mylist.begin(); lcurr != mylist.end(); ++ lcurr )
	{
		fno =  * lcurr;
		le = left_cell_of_face[ fno ];
		re = right_cell_of_face[ fno ];
		cell2coarsegridcell[le] = cnTCell;
		cell2coarsegridcell[re] = cnTCell ++;
	}

	// some cells may not have been merged, merge them now
	for ( int icell = 0; icell < nTotalCell; ++ icell )
	{
		if(cell2coarsegridcell[icell] == -1)
		{
			cell2coarsegridcell[icell] = cnTCell ++;
		}		
	}

	delete [] ratio;
	delete [] unMerged;
	delete [] fmark;

	coarseGrid = new UnstructGrid();
	coarseGrid->InitGrid(new GridID(*fineGrid->GetGridID()), fineGridLevel + 1, fineGrid->GetDim(), UNSTRUCTGRID);

	coarseGrid->SetNTotalCell( cnTCell );
	int nTotalNode = fineGrid->GetNTotalNode();
	coarseGrid->SetNTotalNode( nTotalNode );
	RDouble *x, *y, *z;
	RDouble * xx = fineGrid->GetX();
	RDouble * yy = fineGrid->GetY();
	RDouble * zz = fineGrid->GetZ();

	x = new RDouble[ nTotalNode ];
	y = new RDouble[ nTotalNode ];
	z = new RDouble[ nTotalNode ];

	for ( int inode = 0; inode < nTotalNode; ++ inode )
	{
		x[inode] = xx[inode];
		y[inode] = yy[inode];
		z[inode] = zz[inode];
	}
	coarseGrid->SetX(x);
	coarseGrid->SetY(y);
	coarseGrid->SetZ(z);
}


void Mesh_Agglomeration::AspectRatioModifyBySkewness(RDouble *ratio, RDouble * xcc, RDouble * ycc, RDouble * zcc, const RDouble & limitAngle)
{
	const RDouble PI = 3.141592653589793;

	int nTotalFace = fineGrid->GetNTotalFace();
	int nBoundFace = fineGrid->GetNBoundFace();
	int nTotalCell = fineGrid->GetNTotalCell();
	int n = nBoundFace + nTotalCell;
	int *left_cell_of_face = fineGrid->GetLeftCellOfFace();
	int *right_cell_of_face = fineGrid->GetRightCellOfFace();

	RDouble * xfc = fineGrid->GetFaceCenterX();
	RDouble * yfc = fineGrid->GetFaceCenterY();
	RDouble * zfc = fineGrid->GetFaceCenterZ();

	RDouble * nxs  = fineGrid->GetFaceNormalX();
	RDouble * nys  = fineGrid->GetFaceNormalY();
	RDouble * nzs  = fineGrid->GetFaceNormalZ();

	int * face_number_of_each_cell = fineGrid->GetFaceNumberOfEachCell();
	int ** cell2face = fineGrid->GetCell2Face();

	int LRCell[2];
	RDouble dx, dy, dz, dist, dot;
	for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
	{
		LRCell[0] = left_cell_of_face[iFace];
		LRCell[1] = right_cell_of_face[iFace];


		RDouble minAngle = 100.0;
		for (int iCell = 0; iCell < 2; ++ iCell)
		{
			int cellID = LRCell[iCell];

			int nFace = face_number_of_each_cell[cellID];
			for (int jFace = 0; jFace < nFace; ++ jFace)
			{
				int faceID = cell2face[cellID][jFace];
				if(faceID == iFace) continue;

				dx = xfc[faceID] - xcc[iFace];
				dy = yfc[faceID] - ycc[iFace];
				dz = zfc[faceID] - zcc[iFace];
				dist  = DISTANCE(dx, dy, dz);
				dot = ( nxs[faceID] * dx + nys[faceID] * dy + nzs[faceID] * dz ) / ( dist + TINY );
				if ( dot >  1.0 ) dot =   1.0;
				if ( dot < -1.0 ) dot = - 1.0;

				if(right_cell_of_face[faceID] == cellID)
				{
					dot *= -1;
				}

				RDouble angle  = asin(dot) * 180.0 / PI;
				minAngle = MIN(angle, minAngle);
			}
		}

		RDouble factor = minAngle / 90.0;

		factor = 1.0 - factor;

		ratio[iFace] *= factor;

		if(minAngle < limitAngle) ratio[iFace] = -1;
	}
}

}