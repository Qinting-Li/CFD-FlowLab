inline int  StructGrid::GetMultigridStepI() const 
{ 
	return multigridInfo->GetMultigridStepI(); 
}

inline int  StructGrid::GetMultigridStepJ() const 
{ 
	return multigridInfo->GetMultigridStepJ(); 
}

inline int  StructGrid::GetMultigridStepK() const 
{ 
	return multigridInfo->GetMultigridStepK(); 
}

inline void StructGrid::SetMultigridStepI(int istp) 
{ 
	multigridInfo->SetMultigridStepI(istp);
}

inline void StructGrid::SetMultigridStepJ(int jstp) 
{ 
	multigridInfo->SetMultigridStepJ(jstp);
}

inline void StructGrid::SetMultigridStepK(int kstp) 
{ 
	multigridInfo->SetMultigridStepK(kstp);
}

inline void StructGrid::SetBlankIndex(Int3D * iblank) 
{ 
	this->iblank = iblank; 
}

inline RDouble3D * StructGrid::GetWallDist() const 
{ 
	return walldist; 
}

inline RDouble3D * StructGrid::GetCellCenterX() const 
{ 
	return cellMetrics->GetCellCenterX(); 
}

inline RDouble3D * StructGrid::GetCellCenterY() const 
{ 
	return cellMetrics->GetCellCenterY();
}

inline RDouble3D * StructGrid::GetCellCenterZ() const 
{ 
	return cellMetrics->GetCellCenterZ();
}

inline RDouble4D * StructGrid::GetFaceNormalX() const 
{ 
	return faceMetrics->GetFaceNormalX();
}

inline RDouble4D * StructGrid::GetFaceNormalY() const 
{ 
	return faceMetrics->GetFaceNormalY(); 
}

inline RDouble4D * StructGrid::GetFaceNormalZ() const 
{ 
	return faceMetrics->GetFaceNormalZ();
}

inline RDouble4D * StructGrid::GetFaceArea()    const 
{  
	return faceMetrics->GetFaceArea();
}

inline RDouble4D * StructGrid::GetFaceVectorX() const
{
    return faceMetrics->GetFaceVectorX();
}

inline RDouble4D * StructGrid::GetFaceVectorY() const
{
    return faceMetrics->GetFaceVectorY();
}

inline RDouble4D * StructGrid::GetFaceVectorZ() const
{
    return faceMetrics->GetFaceVectorZ();
}

inline RDouble3D * StructGrid::GetCellVolume()  const 
{ 
	return cellMetrics->GetCellVolume();  
}

inline RDouble3D * StructGrid::GetCellLengthX() const
{
	return cellMetrics->GetCellLengthX();
}

inline RDouble3D * StructGrid::GetCellLengthY() const
{
	return cellMetrics->GetCellLengthY();
}

inline RDouble3D * StructGrid::GetCellLengthZ() const
{
	return cellMetrics->GetCellLengthZ();
}

inline RDouble3D * StructGrid::GetCellVolumeOld()  const 
{ 
	return dynamicGridMetrics->GetCellVolumeOld();
}

//! add by Min Yaobing
inline RDouble5D * StructGrid::GetFaceVector_FD()  const 
{ 
	return faceMetrics->GetFaceVector_FD();  
}

inline RDouble5D * StructGrid::GetFaceNormal_FD()  const 
{ 
	return faceMetrics->GetFaceNormal_FD();  
}

inline RDouble3D * StructGrid::GetCellJacobian()  const 
{ 
	return cellMetrics->GetCellJacobian();  
}
//! add by Min Yaobing

inline RDouble4D * StructGrid::GetFaceNormalVelocity() const 
{ 
	return dynamicGridMetrics->GetFaceNormalVelocity(); 
}

inline void StructGrid::SetNI(int ni) 
{ 
	nodeTopology->SetNI( ni );
}

inline void StructGrid::SetNJ(int nj) 
{ 
	nodeTopology->SetNJ( nj );
}

inline void StructGrid::SetNK(int nk) 
{ 
	nodeTopology->SetNK( nk );
}

inline int  StructGrid::GetNI() const 
{ 
	return nodeTopology->GetNI();
}

inline int  StructGrid::GetNJ() const 
{ 
	return nodeTopology->GetNJ();
}

inline int  StructGrid::GetNK() const 
{ 
	return nodeTopology->GetNK();
}

inline void StructGrid::GetCellIterationIndex(int & ist, int & ied, int & jst, int & jed, int & kst, int & ked, int layer)
{
	int ni = nodeTopology->GetNI();
	int nj = nodeTopology->GetNJ();
	int nk = nodeTopology->GetNK();

	ist = 1 - layer;
	ied = ni - 1 + layer;
	jst = 1 - layer;
	jed = nj - 1 + layer;

	if (nk == 1)
	{
		kst = 1;
		ked = 1;
	}
	else
	{
		kst = 1 - layer;
		ked = nk -1 + layer;
	}
}

inline void StructGrid::GetNodeIterationIndex(int & ist, int & ied, int & jst, int & jed, int & kst, int & ked)
{
	int ni = nodeTopology->GetNI();
	int nj = nodeTopology->GetNJ();
	int nk = nodeTopology->GetNK();

	ist = 1;
	ied = ni;
	jst = 1;
	jed = nj;

	if (nk == 1)
	{
		kst = 1;
		ked = 1;
	}
	else
	{
		kst = 1;
		ked = nk;
	}
}

inline void StructGrid::GetFaceIterationIndex(int & ist, int & ied, int & jst, int & jed, int & kst, int & ked, int iSurface)
{
	int ni = nodeTopology->GetNI();
	int nj = nodeTopology->GetNJ();
	int nk = nodeTopology->GetNK();

	if (iSurface == 1)
	{
		ist = 0;
	    ied = ni-1;
	    jst = 1;
	    jed = nj-1;

	    if (nk == 1)
	    {
	    	kst = 1;
	    	ked = 1;
	    }
	    else
	    {
	    	kst = 1;
	    	ked = nk-1;
	    }
	}
	else if (iSurface == 2)
	{
		ist = 1;
	    ied = ni-1;
	    jst = 0;
	    jed = nj-1;

	    if (nk == 1)
	    {
	    	kst = 1;
	    	ked = 1;
	    }
	    else
	    {
	    	kst = 1;
	    	ked = nk-1;
	    }
	}
	else
	{
		ist = 1;
	    ied = ni-1;
	    jst = 1;
	    jed = nj-1;
	    kst = 0;
	    ked = nk-1;
	}
}

inline void StructGrid::GetND( int & si, int & sj, int & sk ) const
{
	si = nodeTopology->GetNI();
	sj = nodeTopology->GetNJ();
	sk = nodeTopology->GetNK();
}

inline void StructGrid::GetLocalCenter( int gp, int & ip, int & jp, int & kp ) const
{
	int ni = nodeTopology->GetNI();
	int nj = nodeTopology->GetNJ();
	kp = gp / ( ( ni - 1 ) * ( nj - 1 ) );

	int other = gp % ( ( ni - 1 ) * ( nj - 1 ) );

	jp = other / ( ni - 1 );

	ip = other % ( ni - 1 );

	return;
}

inline RDouble3D * StructGrid::GetStructX() const 
{ 
	return nodeTopology->GetStructX();
}

inline RDouble3D * StructGrid::GetStructY() const 
{  
	return nodeTopology->GetStructY();
}

inline RDouble3D * StructGrid::GetStructZ() const 
{  
	return nodeTopology->GetStructZ();
}

inline CompositeBCRegion * StructGrid::GetCompositeBCRegion() const 
{ 
	return composite_bcregion; 
}

inline RDouble4D & StructGrid::GetFaceVelocityX() const 
{
	return * dynamicGridMetrics->GetFaceVelocityX();
}

inline RDouble4D & StructGrid::GetFaceVelocityY() const 
{
	return * dynamicGridMetrics->GetFaceVelocityY();
}

inline RDouble4D & StructGrid::GetFaceVelocityZ() const 
{
	return * dynamicGridMetrics->GetFaceVelocityZ();
}

inline int StructGrid::GetZoneStartPointLabel() const
{
	return oversetGridTopology->GetZoneStartPointLabel();
}

inline int * StructGrid::GetCellTypeContainer() const
{
	return oversetGridTopology->GetCellTypeContainer();
}

inline void StructGrid::SetZoneStartPointLabel( int zoneStartPointLabel )
{
	oversetGridTopology->SetZoneStartPointLabel(zoneStartPointLabel);
}

inline int * StructGrid::GetHingedPointContainer() const
{
	return oversetGridTopology->GetHingedPointContainer();
}

inline RDouble * StructGrid::GetXCoreContainer() const
{
	return oversetGridTopology->GetXCoreContainer();
}

inline RDouble * StructGrid::GetYCoreContainer() const
{
	return oversetGridTopology->GetYCoreContainer();
}

inline RDouble * StructGrid::GetZCoreContainer() const
{
	return oversetGridTopology->GetZCoreContainer();
}

inline int StructGrid::GetNumberOfCores() const
{
	return oversetGridTopology->GetNumberOfCores();
}

inline int * StructGrid::GetILinkPointLabel() const
{
	return oversetGridTopology->GetILinkPointLabel();
}

inline int * StructGrid::GetJLinkPointLabel() const
{
	return oversetGridTopology->GetJLinkPointLabel();
}

inline int * StructGrid::GetKLinkPointLabel() const
{
	return oversetGridTopology->GetKLinkPointLabel();
}

inline void StructGrid::SetZoneStartCenterLabel( int zoneStartCenterLabel )
{
	oversetGridTopology->SetZoneStartCenterLabel(zoneStartCenterLabel);
}