inline void BCRegion::SetBCType(int bcType) 
{
	this->bctype = bcType; 
}

inline int  BCRegion::GetBCType() const 
{
	return bctype;
}

inline int  BCRegion::GetRegionBlock() const 
{ 
	return nbs;
}

inline void BCRegion::SetRegionBlock(int currentSourceZoneID)
{ 
	this->nbs = currentSourceZoneID; 
}

inline const string & BCRegion::GetBCName() const 
{
	return this->bcName;
}

inline void BCRegion::SetBCName(const string & bcName)
{
	this->bcName = bcName;
}

inline void BCRegion::SetFaceMatchingTargetDirIndex( int * dir3dIn )
{
    for( int i = 0; i<3; i++ )
    {
        s_dir3d[i] = dir3dIn[i];
    }
}

inline int * BCRegion::GetStam() 
{ 
	return stam; 
}

inline int * BCRegion::GetRate() 
{ 
	return rate; 
}

inline int * BCRegion::GetFaceDirectionIndex()
{
	return s_lr3d; 
}

inline int * BCRegion::GetFaceMatchingTargetDirIndex()
{
	return s_dir3d; 
}

inline int * BCRegion::GetFaceDirectionIndexOfTargetBlock() 
{
	return t_lr3d;
}

inline void BCRegion::SetIJKRegion(int ist, int ied, int jst, int jed, int kst, int ked)
{
	s_st[0] = ist;
	s_ed[0] = ied;
	s_st[1] = jst;
	s_ed[1] = jed;
	s_st[2] = kst;
	s_ed[2] = ked;
}

inline void BCRegion::GetIJKRegion(int &ist, int &ied, int &jst, int &jed, int &kst, int &ked)
{
	ist    = s_st[0];
	ied    = s_ed[0];
	jst    = s_st[1];
	jed    = s_ed[1];
	kst    = s_st[2];
	ked    = s_ed[2];
}

inline void BCRegion::SetTargetIJKRegion(int ist, int ied, int jst, int jed, int kst, int ked)
{
	t_st[0] = ist;
	t_ed[0] = ied;
	t_st[1] = jst;
	t_ed[1] = jed;
	t_st[2] = kst;
	t_ed[2] = ked;
}

inline int * BCRegion::GetStartPoint()
{
	return s_st; 
}

inline int BCRegion::GetStartPoint(int m) 
{ 
	return s_st[m];
}

inline int BCRegion::GetTargetStart( int m )
{
	return t_st[m];
}

inline int * BCRegion::GetTargetStart()
{
	return t_st; 
}

inline int BCRegion::GetTargetEnd( int m )
{
	return t_ed[m]; 
}

inline int * BCRegion::GetTargetEnd() 
{
	return t_ed; 
}

inline int * BCRegion::GetEndPoint( void ) 
{
	return s_ed;
}

inline int BCRegion::GetEndPoint(int m) 
{
	return s_ed[m]; 
}

inline int  BCRegion::GetTargetRegionBlock() const 
{
	return nbt;
}

inline void BCRegion::SetTargetRegionBlock(int targetZoneID)
{
	this->nbt = targetZoneID;
}

inline int  BCRegion::GetFaceLeftOrRightIndex() const 
{
	return s_lr; 
}

inline void BCRegion::SetFaceLeftOrRightIndex(int s_lr) 
{
	this->s_lr = s_lr;
}

inline int  BCRegion::GetFaceLeftOrRightIndexOfTargetBlock() const 
{
	return t_lr;
}

inline void BCRegion::SetFaceLeftOrRightIndexOfTargetBlock(int t_lr)
{
	this->t_lr = t_lr; 
}

inline int  BCRegion::GetFaceDirection() const 
{
	return s_nd;
}

inline void BCRegion::SetFaceDirection(int s_nd) 
{
	this->s_nd = s_nd;
}

inline int  BCRegion::GetFaceDirectionOfTargetBlock() const 
{
	return t_nd; 
}

inline void BCRegion::SetFaceDirectionOfTargetBlock(int t_nd) 
{
	this->t_nd = t_nd;
}

inline void BCRegion::GetInsideCellIndex(int i, int j, int k, int &is, int &js, int &ks, int nss)
{
	is = i - s_lr3d[0] * (nss-1);
    js = j - s_lr3d[1] * (nss-1);
    ks = k - s_lr3d[2] * (nss-1);
}

inline void BCRegion::GetGhostCellIndex(int i, int j, int k, int &it, int &jt, int &kt, int ntt)
{
	it = i + s_lr3d[0] * ntt;
    jt = j + s_lr3d[1] * ntt;
    kt = k + s_lr3d[2] * ntt;
}

inline void BCRegion::GetBoundaryFaceIndex(int i, int j, int k, int &ib, int &jb, int &kb)
{
	ib = i;
	jb = j;
	kb = k;
	
	if ( s_lr3d[0] == 1 )
	{
		ib = i + 1;
	}
	else if ( s_lr3d[1] == 1 )
	{
		jb = j + 1;
	}
	else if ( s_lr3d[2] == 1 )
	{
		kb = k + 1;
	}
}

inline void CompositeBCRegion::SetNBoundFace(int nBoundFace)
{
	this->nBoundFace = nBoundFace; 
}

inline void CompositeBCRegion::SetNIFace(int nIFace) 
{
	this->nIFace = nIFace; 
}

inline int  CompositeBCRegion::GetNBoundFace() const
{ 
	return nBoundFace;
}

inline int  CompositeBCRegion::GetNIFace() const 
{
	return nIFace;
}

inline void CompositeBCRegion::CreateBCRegion(uint_t nBCRegion)
{
	bcregions = new vector<BCRegion *>(nBCRegion);
}

inline void CompositeBCRegion::SetBCRegion(uint_t iBCRegion, BCRegion * bcregion)
{
	(*bcregions)[iBCRegion]  = bcregion;
}

inline int CompositeBCRegion::GetnBCRegion() const 
{
	return static_cast<int>( bcregions->size()); 
}