inline void Geo_DynamicGridMetrics_Unstruct::SetNodeVelocityX(RDouble * xtn)
{
	this->xtn = xtn;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetNodeVelocityY(RDouble * ytn)
{
	this->ytn = ytn;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetNodeVelocityZ(RDouble * ztn)
{
	this->ztn = ztn;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetFaceNormalVelocity(RDouble * vgn)
{
	this->vgn = vgn;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetCellVolumeOld(RDouble * voln)
{
	this->voln = voln;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetNodeVelocityX() const 
{ 
	return xtn; 
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetNodeVelocityY() const 
{
	return ytn; 
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetNodeVelocityZ() const 
{
	return ztn;  
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetFaceNormalVelocity() const 
{
	return vgn; 
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetCellVolumeOld() const
{
	return voln;
}