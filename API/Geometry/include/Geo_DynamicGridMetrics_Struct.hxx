inline void Geo_DynamicGridMetrics_Struct::DeallocateAll()
{
    if (voln != NULL) { delete voln; voln = NULL; }
    if (xfv  != NULL) { delete xfv;  xfv  = NULL; }
    if (yfv  != NULL) { delete yfv;  yfv  = NULL; }
    if (zfv  != NULL) { delete zfv;  zfv  = NULL; }
    if (vgn  != NULL) { delete vgn;  vgn  = NULL; }
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
    return xfv;
}

inline RDouble4D * Geo_DynamicGridMetrics_Struct::GetFaceVelocityY() const
{
    return yfv;
}

inline RDouble4D * Geo_DynamicGridMetrics_Struct::GetFaceVelocityZ() const
{
    return zfv;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceVelocityX(RDouble4D *xfv)
{
    if (this->xfv && this->xfv != xfv)
    {
        delete this->xfv;
    }
    this->xfv = xfv;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceVelocityY(RDouble4D *yfv)
{
    if (this->yfv && this->yfv != yfv)
    {
        delete this->yfv;
    }
    this->yfv = yfv;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceVelocityZ(RDouble4D *zfv)
{
    if (this->zfv && this->zfv != zfv)
    {
        delete this->zfv;
    }
    this->zfv = zfv;
}

inline void Geo_DynamicGridMetrics_Struct::SetCellVolumeOld(RDouble3D *voln)
{
    if (this->voln && this->voln != voln)
    {
        delete this->voln;
    }
    this->voln = voln;
}

inline void Geo_DynamicGridMetrics_Struct::SetFaceNormalVelocity(RDouble4D *vgn)
{
    if (this->vgn && this->vgn != vgn)
    {
        delete this->vgn;
    }
    this->vgn = vgn;
}