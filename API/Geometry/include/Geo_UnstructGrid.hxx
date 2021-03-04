inline void UnstructGrid::SetLeftCellOfFace(int *left_cell_of_face_in) 
{ 
	this->faceTopology->SetLeftCellOfFace(left_cell_of_face_in);
}

inline void UnstructGrid::SetRightCellOfFace(int *right_cell_of_face_in) 
{
	this->faceTopology->SetRightCellOfFace(right_cell_of_face_in);
}

inline void UnstructGrid::SetNodeNumberOfEachFace(int *node_number_of_each_face)
{
	this->faceTopology->SetNodeNumberOfEachFace(node_number_of_each_face);
}
inline void UnstructGrid::SetFace2Node(int *face2node)
{
	this->faceTopology->SetFace2Node(face2node);
}

inline void UnstructGrid::SetFace2NodeSubscript(int *face2nodeSubscript)
{
    this->faceTopology->SetFace2NodeSubscript(face2nodeSubscript);
}

inline void UnstructGrid::SetCell2CoarseGridCell(int *cell2coarsegridcell) 
{ 
	this->cell2coarsegridcell = cell2coarsegridcell;
}

inline void UnstructGrid::Setfc2cL(int *fc2cL) 
{ 
	this->fc2cL = fc2cL; 
}

inline void UnstructGrid::Setfc2cR(int *fc2cR) 
{ 
	this->fc2cR = fc2cR;
}

inline void UnstructGrid::SetBCRecord(BCRecord **bcr) 
{
	this->bcr = bcr; 
}

inline void UnstructGrid::SetBCRegionUnstruct(BCRegionUnstruct **bcRegionUnStruct) 
{
	this->bcRegionUnstruct = bcRegionUnstruct; 
}

inline void UnstructGrid::SetNBCRegionUnstruct(int nBCRegionUnstruct) 
{
	this->nBCRegionUnstruct = nBCRegionUnstruct; 
}

inline void UnstructGrid::SetBlankIndex(int *iblank) 
{
	this->iblank = iblank; 
}

inline void UnstructGrid::SetProbesNearestCellID(vector <int> probesNearestCellID)
{
    this->probesNearestCellID = probesNearestCellID; 
}

inline int * UnstructGrid::GetLeftCellOfFace() const 
{ 
	return this->faceTopology->GetLeftCellOfFace();
}

inline int * UnstructGrid::GetRightCellOfFace() const
{ 
	return this->faceTopology->GetRightCellOfFace();
}

inline int * UnstructGrid::GetNodeNumberOfEachFace() const
{
	return this->faceTopology->GetNodeNumberOfEachFace();
}

inline int * UnstructGrid::GetFace2Node() const
{
	return this->faceTopology->GetFace2Node();
}

inline RDouble * UnstructGrid::GetFaceCenterX() const
{
	return this->faceMetrics->GetFaceCenterX();
}

inline RDouble * UnstructGrid::GetFaceCenterY() const
{
	return this->faceMetrics->GetFaceCenterY();
}

inline RDouble * UnstructGrid::GetFaceCenterZ() const
{
	return this->faceMetrics->GetFaceCenterZ();
}

inline RDouble * UnstructGrid::GetFaceNormalX() const
{
	return this->faceMetrics->GetFaceNormalX();
}

inline RDouble * UnstructGrid::GetFaceNormalY() const
{
	return this->faceMetrics->GetFaceNormalY();
}

inline RDouble * UnstructGrid::GetFaceNormalZ() const
{
	return this->faceMetrics->GetFaceNormalZ();
}

inline RDouble * UnstructGrid::GetFaceArea() const
{
	return this->faceMetrics->GetFaceArea();
}

inline RDouble * UnstructGrid::GetCellCenterX() const
{
	return cellMetrics->GetCellCenterX();
}

inline RDouble * UnstructGrid::GetCellCenterY() const 
{
	return cellMetrics->GetCellCenterY();
}

inline RDouble * UnstructGrid::GetCellCenterZ() const 
{
	return cellMetrics->GetCellCenterZ();
}

inline RDouble * UnstructGrid::GetCellVolume()  const
{
	return cellMetrics->GetCellVolume();
}

inline RDouble * UnstructGrid::GetNodeVelocityX() const 
{ 
	return dynamicGridMetrics->GetNodeVelocityX();
}

inline RDouble * UnstructGrid::GetNodeVelocityY() const 
{
	return dynamicGridMetrics->GetNodeVelocityY();
}

inline RDouble * UnstructGrid::GetNodeVelocityZ() const 
{
	return dynamicGridMetrics->GetNodeVelocityZ();
}

inline RDouble * UnstructGrid::GetFaceNormalVelocity() const 
{
	return dynamicGridMetrics->GetFaceNormalVelocity();
}

inline RDouble * UnstructGrid::GetCellVolumeOld()  const
{
	return dynamicGridMetrics->GetCellVolumeOld();
}


inline RDouble * UnstructGrid::GetLeastSquareIWT() const
{
	return leastSquareWeights->GetLeastSquareIWT();
}

inline RDouble * UnstructGrid::GetLeastSquareIXX() const
{
	return leastSquareWeights->GetLeastSquareIXX();
}

inline RDouble * UnstructGrid::GetLeastSquareIYY() const
{
	return leastSquareWeights->GetLeastSquareIYY();
}

inline RDouble * UnstructGrid::GetLeastSquareIZZ() const
{
	return leastSquareWeights->GetLeastSquareIZZ();
}

inline RDouble * UnstructGrid::GetLeastSquareIXY() const
{
	return leastSquareWeights->GetLeastSquareIXY();
}

inline RDouble * UnstructGrid::GetLeastSquareIXZ() const
{
	return leastSquareWeights->GetLeastSquareIXZ();
}

inline RDouble * UnstructGrid::GetLeastSquareIYZ() const
{
	return leastSquareWeights->GetLeastSquareIYZ();
}

inline void UnstructGrid::SetLeastSquareIWT(RDouble *iwt)
{
	leastSquareWeights->SetLeastSquareIWT(iwt);
}

inline void UnstructGrid::SetLeastSquareIXX(RDouble *ixx)
{
	leastSquareWeights->SetLeastSquareIXX(ixx);
}

inline void UnstructGrid::SetLeastSquareIYY(RDouble *iyy)
{
	leastSquareWeights->SetLeastSquareIYY(iyy);
}

inline void UnstructGrid::SetLeastSquareIZZ(RDouble *izz)
{
	leastSquareWeights->SetLeastSquareIZZ(izz);
}

inline void UnstructGrid::SetLeastSquareIXY(RDouble *ixy)
{
	leastSquareWeights->SetLeastSquareIXY(ixy);
}

inline void UnstructGrid::SetLeastSquareIXZ(RDouble *ixz)
{
	leastSquareWeights->SetLeastSquareIXZ(ixz);
}

inline void UnstructGrid::SetLeastSquareIYZ(RDouble *iyz)
{
	leastSquareWeights->SetLeastSquareIYZ(iyz);
}

inline void UnstructGrid::SetFaceMark(char *fMark)
{
	leastSquareWeights->SetFaceMark(fMark);
}

inline char * UnstructGrid::GetFaceMark() const
{
	return leastSquareWeights->GetFaceMark();
}

inline int * UnstructGrid::GetCell2CoarseGridCell(void) const 
{
	return cell2coarsegridcell;
}

inline int * UnstructGrid::Getfc2cL(void) const 
{ 
	return fc2cL;
}

inline int * UnstructGrid::Getfc2cR(void) const
{
	return fc2cR;
}

inline BCRecord ** UnstructGrid::GetBCRecord(void) const 
{
	return bcr; 
}

inline int UnstructGrid::GetNBCRegionUnstruct(void) const 
{
	return nBCRegionUnstruct; 
}

inline int * UnstructGrid::GetBlankIndex(void) const 
{
	return iblank; 
}

inline vector <int> UnstructGrid::GetProbesNearestCellID() const 
{
    return probesNearestCellID; 
}

inline RDouble * UnstructGrid::GetWallDist() const
{
	return walldist;
}

inline RDouble * UnstructGrid::GetLamdax() const 
{
	return lamdax;
}

inline RDouble * UnstructGrid::GetLamday() const
{
	return lamday;
}

inline RDouble * UnstructGrid::GetLamdaz() const
{
	return lamdaz;
}

inline int * UnstructGrid::Getknode() const 
{
	return knode;
} 

//! Set the cell to node: node index of each cell.
inline void UnstructGrid::SetCell2Node(int * cell2node)
{
	cellTopology->SetCell2Node(cell2node);
}

inline map<int,int> & UnstructGrid::GetBoundaryPointLabel()
{
    return boundaryPointLabel;
}

//! Set the number of nodes in each cess.
inline void UnstructGrid::SetNodeNumberOfEachCell(int * node_number_of_each_cell)
{
	cellTopology->SetNodeNumberOfEachCell(node_number_of_each_cell);
}

inline void UnstructGrid::SetOversetConfig( OversetConfig * oversetConfigIn )
{
	this->oversetConfig = oversetConfigIn;
}

inline void UnstructGrid::SetAverageVolume(RDouble rIn)
{
	averageVolume = rIn;
}


