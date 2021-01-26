#include "BC_Struct.h"
#include "Glb_Dimension.h"
#include "Math_BasisFunction.h"
#include "OversetGridFactory.h"
#include "BC.h"
namespace PHSPACE
{
BCRegion::BCRegion(int zoneID, uint_t regionID)
{ 
	this->nbs = zoneID; 
	this->nr = regionID; 
    this->bcData =NULL;

	for(int i = 0; i<3; i++)
	{
		t_st[i] = 0;
		t_ed[i] = 0;
		stam[i] = 0;
		rate[i] = 0;
		s_lr3d[i] = 0;
		t_lr3d[i] = 0;
		s_st[i] = 0;
		s_ed[i] = 0;
        s_dir3d[i] = 0;
	}
	s_nd = -1;
	t_nd = -1;
	//s_nd = 0;
	//t_nd = 0;
	s_lr = 0;
	t_lr = 0;
	nbs = 0;
	nbt = 0;
    bctype = 0;
}

BCRegion::~BCRegion()
{
    if (this->bcData != NULL)
    {
        delete bcData;
    }
}

void BCRegion::ComputeRelativeParameters()
{
	vector< int > slam( 3 ), ltam( 3 ), cs( 3 ), ct( 3 );

	for ( int m = 0; m < 3; ++ m )
	{
		if ( ( GetFaceDirection() == -1 && s_st[m] == s_ed[m] ) || GetFaceDirection() == m )
		//if ( s_st[m] == s_ed[m] )
		{
			slam[m] = 0; 

			if ( s_st[m] == 1 )
			{
				cs[m] = 1; 
			}
			else
			{
				cs[m] = -1; 
			}
		}
		else if ( s_st[m] < 0 )
		{
			slam[m] = 1; 

			if ( s_st[m] < s_ed[m] )
			{
				cs[m] = -1; 
			}
			else
			{
				cs[m] = 1; 
			}
		}
		else
		{
			slam[m] = 2; 

			if ( s_st[m] > s_ed[m] )
			{
				cs[m] = -1; 
			}
			else
			{
				cs[m] = 1; 
			}
		}

		if ( ( GetFaceDirectionOfTargetBlock() == -1 && t_st[m] == t_ed[m] ) || GetFaceDirectionOfTargetBlock() == m )
		//if ( t_st[m] == t_ed[m] )
		{
			ltam[0] = m; 

			if ( t_st[m] == 1 )
			{
				ct[m] = -1; 
			}
			else
			{
				ct[m] = 1; 
			}
		}
		else if ( t_st[m] < 0 )
		{
			ltam[1] = m; 

			if ( t_st[m] < t_ed[m] )
			{
				ct[m] = -1; 
			}
			else
			{
				ct[m] = 1; 
			}
		}
		else
		{
			ltam[2] = m; 

			if ( t_st[m] > t_ed[m] )
			{
				ct[m] = -1; 
			}
			else
			{
				ct[m] = 1; 
			}
		}
	}

	for ( int m = 0; m < 3; ++ m )
	{
		int i = slam[m];
		int j = ltam[i];

		stam[m] = j;
		
		if ( cs[m] == ct[j] )
		{
			rate[m] = 1;
		}
		else
		{
			rate[m] = -1;
		}
	}
	return;
}

void BCRegion::ProbeExternalOverlapCells( int * cellTypeContainer, int nci, int ncj, int nck, int & counter )
{
	vector< int > mst( 3 ), med( 3 );
	for ( int m = 0; m < 3; ++ m )
	{
		if ( s_st[m] == s_ed[m] )
		{
			if ( s_st[m] == 1 )
			{
				mst[m] = s_st[m] - 1;
				med[m] = s_st[m];
			}
			else
			{
				mst[m] = s_ed[m] - 2;
				med[m] = s_ed[m] - 1;
			}
		}
		else
		{
			mst[m] = s_st[m] - 1;
			med[m] = s_ed[m] - 1;
		}
	}

	for ( int kCell = mst[2]; kCell <= med[2]; ++ kCell )
	{
		for ( int jCell = mst[1]; jCell <= med[1]; ++ jCell )
		{
			for ( int iCell = mst[0]; iCell <= med[0]; ++ iCell )
			{
				int mm = nci * ncj * kCell + nci * jCell + iCell;

				if( cellTypeContainer[mm] == GENERIC_COLOR )
				{
					cellTypeContainer[mm] = OVERSET_COLOR;
					counter += 1;
				}
			}
		}
	}

	return;
}

BCRegion & BCRegion::operator = (const BCRegion & rhs)
{
    if ( this == &rhs ) return *this;

	this->s_st[0]   = rhs.s_st[0];
	this->s_ed[0]   = rhs.s_ed[0];
	this->s_st[1]   = rhs.s_st[1];
	this->s_ed[1]   = rhs.s_ed[1];
	this->s_st[2]   = rhs.s_st[2];
	this->s_ed[2]   = rhs.s_ed[2];
	this->SetBCType( rhs.GetBCType() );
	this->SetFaceLeftOrRightIndex( rhs.GetFaceLeftOrRightIndex() );
	this->SetFaceDirection( rhs.GetFaceDirection() );
	this->s_lr3d[0] = rhs.s_lr3d[0];
	this->s_lr3d[1] = rhs.s_lr3d[1];
	this->s_lr3d[2] = rhs.s_lr3d[2];

	this->t_st[0]   = rhs.t_st[0];
	this->t_ed[0]   = rhs.t_ed[0];
	this->t_st[1]   = rhs.t_st[1];
	this->t_ed[1]   = rhs.t_ed[1];
	this->t_st[2]   = rhs.t_st[2];
	this->t_ed[2]   = rhs.t_ed[2];

    this->SetRegionBlock( rhs.GetRegionBlock() );
	this->SetTargetRegionBlock( rhs.GetTargetRegionBlock() );
	this->SetFaceDirectionOfTargetBlock ( rhs.GetFaceDirectionOfTargetBlock() );
	this->SetFaceLeftOrRightIndexOfTargetBlock(  rhs.GetFaceLeftOrRightIndexOfTargetBlock() );
	this->t_lr3d[0] = rhs.t_lr3d[0];
	this->t_lr3d[1] = rhs.t_lr3d[1];
    this->t_lr3d[2] = rhs.t_lr3d[2];
	this->stam[0] = rhs.stam[0];
	this->stam[1] = rhs.stam[1];
    this->stam[2] = rhs.stam[2];
	this->rate[0] = rhs.rate[0];
	this->rate[1] = rhs.rate[1];
    this->rate[2] = rhs.rate[2];
    this->bcName  = rhs.GetBCName();

    return *this;
}

void BCRegion::GetNormalizeIJKRegion(int &ist, int &ied, int &jst, int &jed, int &kst, int &ked)
{
    ist = MIN( ABS(this->s_st[0]), ABS(this->s_ed[0]) );
    ied = MAX( ABS(this->s_st[0]), ABS(this->s_ed[0]) );
    jst = MIN( ABS(this->s_st[1]), ABS(this->s_ed[1]) );
    jed = MAX( ABS(this->s_st[1]), ABS(this->s_ed[1]) );
    kst = MIN( ABS(this->s_st[2]), ABS(this->s_ed[2]) );
    ked = MAX( ABS(this->s_st[2]), ABS(this->s_ed[2]) );
}

void BCRegion::ComputeMultiGridIJKRegion(int imin, int imax, int jmin, int jmax, int kmin, int kmax, int istp, int jstp, int kstp, int nsurf)
{
	s_st[0] = ( imin - 1 ) / istp + 1;
	s_ed[0] = MAX(imax / istp, 1);

	s_st[1] = ( jmin - 1 ) / jstp + 1;
	s_ed[1] = MAX(jmax / jstp, 1);

	s_st[2] = ( kmin - 1 ) / kstp + 1;
	s_ed[2] = MAX(kmax / kstp, 1);

    s_ed[nsurf-1] = s_st[nsurf-1];
}

void BCRegion::InitFaceDirectionIndex()
{
	this->s_lr3d[0] = 0;
	this->s_lr3d[1] = 0;
    this->s_lr3d[2] = 0;
	//! Here do some special handling.
	this->s_lr3d[ s_nd ] = s_lr;
}

void BCRegion::ProcessBCInfo()
{
    using namespace PHSPACE;

	int s_st_local[3],s_ed_local[3];
	int t_st_local[3],t_ed_local[3];
	int s_fix;

    int nDim = GetDim();

	/*for ( int m = 0; m < 3; ++ m )
	{
		s_st[m] = this->s_st[m];
		s_ed[m] = this->s_ed[m];
	}

	for ( int m = 0; m < 3; ++ m )
	{
		this->s_st[m] = MIN( ABS(s_st[m]), ABS(s_ed[m]) );
		this->s_ed[m] = MAX( ABS(s_st[m]), ABS(s_ed[m]) );

		this->s_lr3d[m] = 0;
    }*/


	for ( int m = 0; m < 3; ++ m )
	{
		s_st_local[m] = this->s_st[m];
		s_ed_local[m] = this->s_ed[m];

		t_st_local[m] = this->t_st[m];
		t_ed_local[m] = this->t_ed[m];
	}

    if (IsInterface(bctype) && nDim == THREE_D)
    {
        for (int m = 0; m < 3; ++ m)
        {
            if (s_st_local[m] == s_ed_local[m])
            {
                for (int n = 0; n < 3; ++ n)
                {
                    if (t_st_local[n] == t_ed_local[n])
                    {
                        if (((s_st_local[m] == 1) && (t_st_local[n] == 1)) || ((s_st_local[m] != 1) && (t_st_local[n] != 1)))
                        {
                            s_dir3d[m] = - n - 1;
                        }
                        else
                        {
                            s_dir3d[m] = n + 1;
                        }
                        break;
                    }
                    else
                    {
                        continue;
                    }
                }
            }
            else
            {
                for (int n = 0; n < 3; ++ n)
                {
                    if ((t_st_local[n] == t_ed_local[n]) || (t_st_local[n] * s_st_local[m] < 0))
                    {
                        continue;
                    }
                    else
                    {
                        if (((s_st_local[m] < s_ed_local[m]) && (t_st_local[n] < t_ed_local[n])) || ((s_st_local[m] > s_ed_local[m]) && (t_st_local[n] > t_ed_local[n])))
                        {
                            s_dir3d[m] = n + 1;
                        }
                        else
                        {
                            s_dir3d[m] = - n - 1;
                        }
                        break;
                    }
                }
            }
        }
    }

	for ( int m = 0; m < 3; ++ m )
	{
		if ( ABS(s_st_local[m]) > ABS(s_ed_local[m]) )
		{
			this->s_st[m] = ABS(s_ed_local[m]);
			this->s_ed[m] = ABS(s_st_local[m]);

			this->t_st[m] = ABS(t_ed_local[m]);
			this->t_ed[m] = ABS(t_st_local[m]);
		}
		else
		{
			this->s_st[m] = ABS(s_st_local[m]);
			this->s_ed[m] = ABS(s_ed_local[m]);

			this->t_st[m] = ABS(t_st_local[m]);
			this->t_ed[m] = ABS(t_ed_local[m]);
		}

		this->s_lr3d[m] = 0;
    }


    for ( int m = 0; m < nDim; ++ m )
    {
		if ( ( GetFaceDirection() == -1 && s_st_local[m] == s_ed_local[m] ) || GetFaceDirection() == m )
		//if ( s_st[m] == s_ed[m] )
		{
			this->SetFaceDirection( m );
			if ( s_st_local[m] == 1 )
			{
				this->SetFaceLeftOrRightIndex( - 1 );
				this->s_lr3d[m] = - 1;
			}
			else
			{
				this->SetFaceLeftOrRightIndex( 1 );
				this->s_lr3d[m] = 1;
			}
			s_fix = s_st_local[m];

			if ( this->GetFaceLeftOrRightIndex() == 1 )
			{
				s_st_local[m] -= 1;
				s_ed_local[m] -= 1;
				this->s_st[m] -= 1;
				this->s_ed[m] -= 1;
			}
		}
		else
		{
			this->s_ed[m] -= 1;
		}
	}

	if ( ! IsInterface( bctype ) || nDim == TWO_D ) return;

	int tst[3], ted[3];
	for ( int m = 0; m < 3; ++ m )
	{
		tst[m] = t_st_local[m];
		ted[m] = t_ed_local[m];
	}

	for ( int m = 0; m < 3; ++ m )
	{
		this->t_st[m] = min( abs( tst[m] ), abs( ted[m] ) );
		this->t_ed[m] = max( abs( tst[m] ), abs( ted[m] ) );
		this->t_lr3d[m] = 0;
    }

    for ( int m = 0; m < 3; ++ m )
    {
		if ( ( GetFaceDirectionOfTargetBlock() == -1 && tst[m] == ted[m] ) || GetFaceDirectionOfTargetBlock() == m )
		//if ( tst[m] == ted[m] )
		{
			//this->t_nd = m;
			this->SetFaceDirectionOfTargetBlock( m );
			if ( tst[m] == 1 )
			{
				this->t_lr = -1;
				this->t_lr3d[m] = -1;
			}
			else
			{
				this->t_lr = 1;
				this->t_lr3d[m] = 1;
			}

			if ( this->t_lr == 1 )
			{
				this->t_st[m] -= 1;
				this->t_ed[m] -= 1;
			}
		}
		else
		{
			this->t_ed[m] -= 1;
		}
	}

	return;
}

CompositeBCRegion::CompositeBCRegion(int zoneID)
{
    this->zoneID = zoneID;
	bcregions   = 0;
	id_i  = 0;
	id_j  = 0;
	id_k  = 0;
	face_of_bcregion = 0;

    InitBoundaryName();
}

CompositeBCRegion::~CompositeBCRegion()
{
	if ( bcregions )
	{
		for ( std::size_t i = 0; i < bcregions->size(); ++ i )
		{
			delete (*bcregions)[i];
		}
	}
	delete bcregions;

    DeleteIndex();
}



void CompositeBCRegion::CopyCompositeBCRegion(CompositeBCRegion * rightHandSide)
{
	int nBCRegion = rightHandSide->GetnBCRegion();

	CreateBCRegion(nBCRegion);

	for ( uint_t ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion )
	{
		BCRegion * bcregion = new BCRegion(zoneID, ibcregion);

		this->SetBCRegion(ibcregion, bcregion);
		*bcregion = *rightHandSide->GetBCRegion(ibcregion);
	}
}

BCRegion * CompositeBCRegion::GetBCRegion(uint_t iBCRegion) const 
{
    vector<BoundaryCondition *> * globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    Data_Param * bcData = (*bcregions)[iBCRegion]->GetBCDataParam();
    if (bcData == NULL)
    {
        if (globalBCList != NULL)
        {
            vector<BoundaryCondition *>::iterator iter;
            string bcRegionName = (*bcregions)[iBCRegion]->GetBCName();
            for (iter = globalBCList->begin(); iter != globalBCList->end(); iter++)
            {
                string globalBCName = (*iter)->GetBCName();
                Data_Param * globalBCParam = (*iter)->GetBCParamDataBase();
                if (globalBCName == bcRegionName)
                {
                    (*bcregions)[iBCRegion]->SetBCDataParam(globalBCParam);
                }
            }
        }
        else
        {
            Data_Param * bcData = new Data_Param();
            (*bcregions)[iBCRegion]->SetBCDataParam(bcData);
        }
        return (*bcregions)[iBCRegion];
    }

    return (*bcregions)[iBCRegion]; 
}

void CompositeBCRegion::DeleteIndex()
{
	delete [] id_i;
	delete [] id_j;
	delete [] id_k;
	delete [] face_of_bcregion;

    //! The reason for setting zero is that DeleteIndex may be called many times.
	id_i  = 0;
	id_j  = 0;
	id_k  = 0;
	face_of_bcregion = 0;
}

void CompositeBCRegion::DecodeIJK(int index, int &i, int &j, int &k)
{
	i = id_i[ index ];
	j = id_j[ index ];
	k = id_k[ index ];
}

void CompositeBCRegion::GetSourceIndexIJK(int iface, int ipos, int &i, int &j, int &k)
{
	int ibcregion = face_of_bcregion[iface];
	BCRegion * bcregion = this->GetBCRegion(ibcregion);

	int is,js,ks;
	DecodeIJK(iface, is, js, ks);

	int * s_lr3d = bcregion->GetFaceDirectionIndex();

	i = is - s_lr3d[0] * ( ipos - 1 );
	j = js - s_lr3d[1] * ( ipos - 1 );
	k = ks - s_lr3d[2] * ( ipos - 1 );
}

void CompositeBCRegion::GetSourceIndexIJK_Nsurf_LR(int iface, int ipos, int &i, int &j, int &k, int &nsurf, int &s_lr)
{
	if (!face_of_bcregion) SetIFaceInfo();

	int ibcregion = face_of_bcregion[iface];
	BCRegion * bcregion = this->GetBCRegion(ibcregion);

	int is,js,ks;
	DecodeIJK(iface, is, js, ks);

	int * s_lr3d = bcregion->GetFaceDirectionIndex();

	i = is - s_lr3d[0] * ( ipos - 1 );
	j = js - s_lr3d[1] * ( ipos - 1 );
	k = ks - s_lr3d[2] * ( ipos - 1 );

	s_lr = bcregion->GetFaceLeftOrRightIndex();
	nsurf = bcregion->GetFaceDirection()+1;
}

void CompositeBCRegion::GetTargetIndexIJK(int iface, int ipos, int &i, int &j, int &k)
{
	int ibcregion = face_of_bcregion[iface];
	BCRegion * bcregion = this->GetBCRegion(ibcregion);

	int is,js,ks;
	DecodeIJK(iface, is, js, ks);

	int * s_lr3d = bcregion->GetFaceDirectionIndex();

	i = is + s_lr3d[0] * ( ipos );
	j = js + s_lr3d[1] * ( ipos );
	k = ks + s_lr3d[2] * ( ipos );
}

void CompositeBCRegion::GetSourceIndexIJK_fornodetransfer(int iface, int &id, int &jd, int &kd, int &i_lr, int &j_lr, int &k_lr)
{
	int ibcregion = face_of_bcregion[iface];
	BCRegion * bcregion = this->GetBCRegion(ibcregion);
    
	int * s_lr3d = bcregion->GetFaceDirectionIndex();

    i_lr = s_lr3d[0];
    j_lr = s_lr3d[1];
    k_lr = s_lr3d[2];
    
    id = 0;
	jd = 0;
    kd = 0;

	if ( s_lr3d[0] == 1 )
	{
		id = 1;
	}
	else if ( s_lr3d[1] == 1 )
	{
		jd = 1;
	}
	else if ( s_lr3d[2] == 1 )
	{
		kd = 1;
	}
}

void CompositeBCRegion::GetTargetIndexIJK_fornodetransfer(int iface, int &id, int &jd, int &kd, int &i_lr, int &j_lr, int &k_lr)
{
	int ibcregion = face_of_bcregion[iface];
	BCRegion * bcregion = this->GetBCRegion(ibcregion);

	int * s_lr3d = bcregion->GetFaceDirectionIndex();

    i_lr = s_lr3d[0];
    j_lr = s_lr3d[1];
    k_lr = s_lr3d[2];

    id = 0;
	jd = 0;
    kd = 0;

	if ( s_lr3d[0] == -1 )
	{
		id = 1;
	}
	else if ( s_lr3d[1] == -1 )
	{
		jd = 1;
	}
	else if ( s_lr3d[2] == -1 )
	{
		kd = 1;
	}
}



void CompositeBCRegion::SetIFaceInfo()
{
	using namespace PHSPACE;

    DeleteIndex();

	int nBCRegion = this->GetnBCRegion();
	int imin,imax,jmin,jmax,kmin,kmax,bctype;

	nIFace = 0;
    nBoundFace = 0;
	for ( uint_t ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion )
	{
		BCRegion * bcregion = this->GetBCRegion(ibcregion);
		bcregion->GetIJKRegion(imin,imax,jmin,jmax,kmin,kmax);

        nBoundFace += ( imax - imin + 1 ) * ( jmax - jmin + 1 ) * ( kmax - kmin + 1 );

		bctype = bcregion->GetBCType();
		if ( IsInterface(bctype) )
		{
			nIFace += ( imax - imin + 1 ) * ( jmax - jmin + 1 ) * ( kmax - kmin + 1 );
		}
	}

	if ( nIFace == 0 ) return;

	id_i  = new int [ nIFace ];
	id_j  = new int [ nIFace ];
    id_k  = new int [ nIFace ];
	face_of_bcregion = new int [ nIFace ];

	int iiface = 0;
	for ( int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion )
	{
		BCRegion * bcregion = this->GetBCRegion(ibcregion);
		bcregion->GetIJKRegion(imin,imax,jmin,jmax,kmin,kmax);
		bctype = bcregion->GetBCType();

        if (IsInterface(bctype))
		{
			for ( int k = kmin; k <= kmax; ++ k )
			{
				for ( int j = jmin; j <= jmax; ++ j )
				{
					for ( int i = imin; i <= imax; ++ i )
					{
						id_i [iiface] = i;
						id_j [iiface] = j;
						id_k [iiface] = k;
						face_of_bcregion[iiface] = ibcregion;
						++ iiface;
					}
				}
			}
        }
	}
}

void CompositeBCRegion::ProcessBCInfo()
{
    using namespace PHSPACE;
	int nBCRegion = this->GetnBCRegion();

	for ( uint_t ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion )
	{
        BCRegion *bcregion = this->GetBCRegion(ibcregion);
		bcregion->ProcessBCInfo();
	}
}

void CompositeBCRegion::InitBoundaryName()
{
    string boundaryName[100] = {"NO_BOUNDARY_CONDITION", "EXTRAPOLATION", "SOLID_SURFACE", "SYMMETRY", 
        "FARFIELD", "INFLOW", "OUTFLOW", "POLE", "Wall_8", "Wall_9", "Wall_10", "Wall_11", 
        "Wall_12", "Wall_13", "Wall_14", "Wall_15", "Wall_16", "Wall_17", "Wall_18", "Wall_19",
        "Wall_20", "Wall_21", "Wall_22", "Wall_23", "Wall_24", "Wall_25", "Wall_26", "Wall_27", 
        "Wall_28", "Wall_29", "Wall_30", "Wall_31", "Wall_32", "Wall_33", "Wall_34", "Wall_35",
        "Wall_36", "Wall_37", "Wall_38", "Wall_39", "Wall_40", "Wall_41", "Wall_42", "Wall_43",
        "Wall_44", "Wall_45", "Wall_46", "Wall_47", "Wall_48", "Wall_49", "Wall_50"};

    boundaryName[52] = "PRESSURE_INLET";
    boundaryName[62] = "PRESSURE_OUTLET";
    boundaryName[71] = "POLE1";
    boundaryName[72] = "POLE2";
    boundaryName[73] = "POLE3";
    //! Add patched face boundary condition by Wang Ziwei.
    boundaryName[81] = "PATCH_FACE1";
    boundaryName[82] = "PATCH_FACE2";
    boundaryName[83] = "PATCH_FACE3";
    boundaryName[84] = "PATCH_FACE4";
    boundaryName[85] = "PATCH_FACE5";
    boundaryName[86] = "PATCH_FACE6";
    boundaryName[87] = "PATCH_FACE7";
    boundaryName[88] = "PATCH_FACE8";
    boundaryName[89] = "PATCH_FACE9";
    boundaryName[90] = "PATCH_FACE10";
    boundaryName[99] = "PERIODIC";

    for (size_t iName = 0; iName < boundaryName->size(); ++ iName)
    {
        this->boundaryName[iName] = boundaryName[iName];
    }
}

void GetBCFaceIDX(int *s_lr3d, int &id, int &jd, int &kd)
{
	id = 0;
	jd = 0;
    kd = 0;

	if ( s_lr3d[0] == 1 )
	{
		id = 1;
	}
	else if ( s_lr3d[1] == 1 )
	{
		jd = 1;
	}
	else if ( s_lr3d[2] == 1 )
	{
		kd = 1;
	}
}

}
