#include "GridBase.h"
using namespace std;
using namespace PHSPACE;
namespace Reorder
{
	//Reorder face label for reducing GPU global memory accessing overheads
	//In CGNSBase::ComputeGridConn after CallReOrderCellLabel
	void CallReOrderFaceLabel(Base_Grid_Conn * gconn);
	//Reorder cell label for reducing bandwidth
	//In CGNSBase::ComputeGridConn after PostProcess
	void CallReOrderCellLabel(Base_Grid_Conn * gconn);
	//Evaluate the bandwidth in CGNSBase::ComputeGridConn
	void CallBandwidthEvaluate(Base_Grid_Conn *gconn);
	//Evaluate the bandwidth in GPUGeomInfoAllocCopy of BasicDeviceVariables
	void CallBandwidthFTS(const int nTotalCell, vector<int> * cell2Cell, const int * ngb_cell_number_of_each_cell);
	//Get the current level of vertex rStart
	//Input: rStart, cell_number_of_each_cell, cell2cell, offsetCell2Cell, ReOrderCell
	//Output: ReOrderCell, numSortCur, cellIDSortCur
	void GetCurrentLevel(int rStart, const int * cell_number_of_each_cell, const int * cell2cell, const int * offsetCell2Cell, int * ReOrderCell, int & numSortCur, int * cellIDSortCur);
	//Get the next level of current level
	//Input:
	//Output:
	void GetNextLevel(const int & numSortCur, const int * cellIDSortCur, const int * cell_number_of_each_cell, const int * cell2cell, const int * offsetCell2Cell, int * ReOrderCell, int & numSortNext, int * cellIDSortNext);
	//copy cellIDSortNext into cellIDSortCur	
	void SetCellIDSortCurByCellIDSortNext(const int & numSortNext, const int * cellIDSortNext, int & numSortCur, int * cellIDSortCur);
	//Reset cellIDSortCur or cellIDSortNext by -1
	void ResetCellSort(int & numSort, int * cellIDSort);
	//Output of cellIDSortCur or cellIDSortNext
	void TestCellSort(const int numSort, const int * cellIDSort);
	//Find the smallest degree in cellIDSortCur
	int FindMinCellIDSortCur(const int numSortCur, const int * cellIDSortCur, const int * cell_number_of_each_cell);
	//Get current level of rStart for reOrder, the difference from  GetCurrentLevel	is that cancel the set of ReOrderCell by 1.
	void GetCurrentLevelForReOrder(int rStart, const int * cell_number_of_each_cell, const int * cell2cell, const int * offsetCell2Cell, int * ReOrderCell, int & numSortCur, int * cellIDSortCur);
	//Arrange cellIDSortCur by Ascending order
	void AscendingCellIDSortCurByDegree(const int numSortCur, const int * cell_number_of_each_cell, int * cellIDSortCur);
	//Test degree of cellIDSortCur from 0 to numSortCur
	void TestCellIDSortCurDegree(const int numSortCur, const int * cellIDSortCur, const int * cell_number_of_each_cell);
	//Check cell2cell by 3 rules:
	//rule1 one cell's neighbor cells should not contain itself.
	//rule2 0<= one cell's neighbor cells < nTotalCell
	//rule3 one cell's neighbor cells should be different.
	void CheckCell2Cell(const int nTotalCell, const int * cell_number_of_each_cell, const int * offsetCell2Cell, const int * cell2cell);
	//Reorder cells in cellIDSortCur
	void ReOrderLevelCur(const int numSortCur, const int * cellIDSortCur, int & numVertex, int * ReOrderCell);

	void GetNextLevelForReOrder(const int & numSortCur, const int * cellIDSortCur, const int * cell_number_of_each_cell, const int * cell2cell, const int * offsetCell2Cell, int * ReOrderCell, int & numSortNext, int * cellIDSortNext);
	//check rules for cell index
	//rule1 0 <= NewCell < nTotalCell
	//rule2 cell index should be different
	void CheckNewCell(const int nTotalCell, const int * NewCell);
	//compute face_number_of_each_cell and nTotalCellFaces
	void SetFaceNumberOfEachCell(const int nTotalFace, const int nTotalCell, const int * left_cell_of_face, const int * right_cell_of_face, int * face_number_of_each_cell, int & nTotalCellFaces);
	//set offsetCell2Face by face_number_of_each_cell
	void SetOffSetCell2Face(const int nTotalCell, const int * face_number_of_each_cell, int * offsetCell2Face);
	//set cell2face by left_cell_of_face, right_cell_of_face, offsetCell2Face
	//It should be noted that boundary faces are not considered here.
	void SetCell2Face(const int nTotalFace, const int nTotalCell, const int * left_cell_of_face, const int * right_cell_of_face, int * offsetCell2Face, int * cell2face);
	//set cell2Cell
	void SetCell2Cell(const int nTotalCell, const int * face_number_of_each_cell, const int * offsetCell2Face, const int * cell2face, const int * left_cell_of_face, const int * right_cell_of_face, int * cell2cell);
	//compute bandwidth of each cell
	void computeBandwidth(const int nTotalCell, const int * face_number_of_each_cell, const int * offsetCell2Face, const int * cell2cell);
	//Update face_number_of_each_face by swap method
	//Intput: nTotalCell, face_number_of_each_face, reflectCell
	//Output:SwapFace_number_of_each_cell
	//The function is used to provide a swap method to update a cell related varible
	//after cell reorder. In the future, cell2node can be update by the method
	void SwapFaceNumberOfEachFace(const int nTotalCell, const int * face_number_of_each_cell, const int * reflectCell, const int * NewCell, const int * NewFace_number_of_each_cell);
	//check swap method for face_number_of_each_face
	//After swap, OrgCell=reflectCell, SwapCell=NewCell, SwapFace_number_of_each_cell=NewFace_number_of_each_cell
	void CheckSwapFaceNumber(const int nTotalCell, const int * OrgCell, const int * SwapCell, const int * reflectCell, const int * NewCell, const int * SwapFace_number_of_each_cell, const int * NewFace_number_of_each_cell);
	//Reorder cell label from nBoundFace to nTotalFace just with the order of faces in cell2face
	void cellReorderByCell2FaceOrder(const int nTotalCell, const int nTotalFace, const int * offsetCell2Face, const int * face_number_of_each_cell, const int * cell2face, const int nBoundFace, int labelFace, int * NewFaces);
	//a caller for MetricLeftRightCellOfFace
	void CallMetricLeftRightCellOfFace(Base_Grid_Conn * gconn);
	//metric for computing left_cell_of_face and right_cell_of_face for GPU 
	void MetricLeftRightCellOfFace(const int nTotalFace, const int nBoundFace, const int * left_cell_of_face, const int * right_cell_of_face);
}
