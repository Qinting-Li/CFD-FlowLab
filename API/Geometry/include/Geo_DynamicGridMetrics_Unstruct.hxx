inline void Geo_DynamicGridMetrics_Unstruct::SetFaceVelocityX(RDouble *xfv)
{
    this->xfv = xfv;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetFaceVelocityY(RDouble *yfv)
{
    this->yfv = yfv;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetFaceVelocityZ(RDouble *zfv)
{
    this->zfv = zfv;
}
inline void Geo_DynamicGridMetrics_Unstruct::SetFaceNormalVelocity(RDouble *vgn)
{
    this->vgn = vgn;
}

inline void Geo_DynamicGridMetrics_Unstruct::SetCellVolumeOld(RDouble *voln)
{
    this->voln = voln;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetFaceVelocityX() const
{
    return xfv;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetFaceVelocityY() const
{
    return yfv;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetFaceVelocityZ() const
{
    return zfv;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetFaceNormalVelocity() const
{
    return vgn;
}

inline RDouble * Geo_DynamicGridMetrics_Unstruct::GetCellVolumeOld() const
{
    return voln;
}