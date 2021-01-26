inline void Geo_DynamicGridMetrics_Struct::DeallocateAll()
{
	if( voln ) { delete voln; voln = 0; }
	if( xtn ) { delete xtn; xtn = 0; }
	if( ytn ) { delete ytn; ytn = 0; }
	if( ztn ) { delete ztn; ztn = 0; }
	if( vgn ) { delete vgn; vgn = 0; }
}

inline RDouble3D * Geo_DynamicGridMetrics_Struct::GetCellVolumeOld() const
{
	return voln;
}

inline RDouble4D * Geo_DynamicGridMetrics_Struct::GetFaceNormalVelocity() const
{
	return vgn;
}

inline RDouble4D * Geo_DynamicGridMetrics_Struct::GetFaceVelocityX() const
{
	return xtn;
}

inline RDouble4D * Geo_DynamicGridMetrics_Struct::GetFaceVelocityY() const
{
	return ytn;
}

inline RDouble4D * Geo_DynamicGridMetrics_Struct::GetFaceVelocityZ() const
{
	return ztn;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceVelocityX(RDouble4D * xtn)
{
	if(this->xtn && this->xtn != xtn)
	{
		delete this->xtn;
	}
	this->xtn = xtn;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceVelocityY(RDouble4D * ytn)
{
	if(this->ytn && this->ytn != ytn)
	{
		delete this->ytn;
	}
	this->ytn = ytn;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceVelocityZ(RDouble4D * ztn)
{
	if(this->ztn && this->ztn != ztn)
	{
		delete this->ztn;
	}
	this->ztn = ztn;
}

inline void Geo_DynamicGridMetrics_Struct::SetCellVolumeOld(RDouble3D * voln)
{
	if(this->voln && this->voln != voln)
	{
		delete this->voln;
	}
	this->voln = voln;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceNormalVelocity(RDouble4D * vgn)
{
	if(this->vgn && this->vgn != vgn)
	{
		delete this->vgn;
	}
	this->vgn = vgn;
}