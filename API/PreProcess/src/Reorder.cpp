#include "Reorder.h"
namespace Reorder
{
	void CallBandwidthFTS(const int nTotalCell, vector<int> * cell2Cell, const int * ngb_cell_number_of_each_cell){
		printf("Program is running in CallBandwidthEvaluate\n");
		int * bandwidth_of_each_cell = new int[nTotalCell];
		for (int icell = 0; icell < nTotalCell; icell++){
			int numCells = ngb_cell_number_of_each_cell[icell];
			int localMaxBand = 0;
			for (int ngbCell = 0; ngbCell < numCells; ngbCell++){
				int bandDiff = abs(icell - cell2Cell[icell][ngbCell]);
				if (bandDiff > localMaxBand) localMaxBand = bandDiff;
			}
			bandwidth_of_each_cell[icell] = localMaxBand;
		}

		//statistics of bandwidth_of_each_cell
		int maxGlobalBand = 0; //maximum bandwidth
		int maxCellID = 0; //cellID of maximum bandwidth
		int bandDist[33]; //distribution of bandwidth from 1 to 32 and more
		for (int i = 0; i < 33; i++){
			bandDist[i] = 0;
		}

		for (int icell = 0; icell < nTotalCell; icell++){
			int bandCell = bandwidth_of_each_cell[icell];
			if (bandCell <= 0) {
				printf("Error: bandwidth<=0 on cell %d, bandwidth=%d\n", icell, bandCell);
				exit(1);
			}
			if (bandCell > maxGlobalBand) {
				maxGlobalBand = bandCell;
				maxCellID = icell;
			}
			if (bandCell > 32) bandDist[32]++;
			else bandDist[bandCell-1]++;
		}
		
		printf("maxGlobalBand = %d on cellID %d\n", maxGlobalBand, maxCellID);
		int nTotalDist = 0;
		for (int i = 0; i < 33; i++){
			nTotalDist += bandDist[i];
			printf("bandDist %d: %d in nTotalCells\n", i+1, bandDist[i]);
		}

		printf("nTotalDist = %d, nTotalCell = %d\n", nTotalDist, nTotalCell);

	}

	void CallBandwidthEvaluate(Base_Grid_Conn *gconn){
		printf("Program is running in CallBandwidthEvaluate\n");
		cgsize_t nTotalFace = gconn->GetNTotalFace();
		cgsize_t nBoundFace = gconn->GetNBoundFace();
		cgsize_t nTotalCell = gconn->GetNTotalCell();
		printf("nTotalFace = %d, nBoundFace = %d, nTotalCell = %d\n", nTotalFace, nBoundFace, nTotalCell);
		//How to get left_cell_of_face and right_cell_of_face is from FYGridBase.cpp
		vector<cgsize_t> &left_cell_of_face  = *gconn->GetLeftCellOfFace();
    		vector<cgsize_t> &right_cell_of_face = *gconn->GetRightCellOfFace();
		//It is difficult to create cell2face, because the number of faces belonging to a face is not determined. So, I decide to use the familiar schemes
		int * cell2face;

		int * cell2cell;
		//offset of icell's cell faces in cell2face
		int * offsetCell2Face = new int[nTotalCell];
		int nTotalCellFaces = 0;

		int * face_number_of_each_cell = new int[nTotalCell];

		int * bandwidth_of_each_cell = new int[nTotalCell];

		//Initialization of face_number_of_each_cell
		printf("Initialization\n");
		for (int icell = 0; icell < nTotalCell; icell++){
			face_number_of_each_cell[icell] = 0;
			bandwidth_of_each_cell[icell] = 0;
			offsetCell2Face[icell] = 0;
		}
		//set face_number_of_each_cell and nTotalCellFaces
		printf("set face_number_of_each_cell and nTotalCellFaces\n");
		for (int iface = 0; iface < nTotalFace; iface++){
			//for left_cell_of_face
			int cellLabelLeft = left_cell_of_face[iface];
			if ((cellLabelLeft < 0)||(cellLabelLeft >= nTotalCell)) {
				printf("Error: cellLabelLeft = %d\n", cellLabelLeft);
				exit(1);
			}
			//face_number_of_each_cell[cellLabel]++;
			//nTotalCellFaces++;
			//for right_cell_of_face
			int cellLabelRight = right_cell_of_face[iface];
			if (cellLabelRight == -1) continue; //boundy face's right neighbor

			if ((cellLabelRight < 0)||(cellLabelRight >= nTotalCell)) {
				printf("Error: cellLabelRight = %d\n", cellLabelRight);
				exit(1);
			}
			face_number_of_each_cell[cellLabelLeft]++;
			nTotalCellFaces++;
			face_number_of_each_cell[cellLabelRight]++;
			nTotalCellFaces++;
		}		
		printf("set the size of cell2cell and cell2face\n");
		//test begin
			//printf("On cell %d, face_number_of_each_cell[%d]=%d\n", 0, 0, face_number_of_each_cell[0]);
		//test end
		//set the size of cell2cell and cell2face
		cell2face = new int[nTotalCellFaces];
		cell2cell = new int[nTotalCellFaces];
		//set offsetCell2Face
		offsetCell2Face[0] = 0;
		for (int icell = 1; icell < nTotalCell; icell++){
			offsetCell2Face[icell] = offsetCell2Face[icell-1] + face_number_of_each_cell[icell-1];
		}
		//set cell2face
		for (int iface = 0; iface < nTotalFace; iface++){
			//for left_cell_of_face
			int cellLabelLeft = left_cell_of_face[iface];
			if ((cellLabelLeft < 0)||(cellLabelLeft >= nTotalCell)) {
				printf("Error: cellLabelLeft = %d\n", cellLabelLeft);
				exit(1);
			}
			//int cell2faceOffset = offsetCell2Face[cellLabel];
			//cell2face[cell2faceOffset] = iface;
			//offsetCell2Face update for the next set of face
			//offsetCell2Face[cellLabel]++;
			//for right_cell_of_face
			int cellLabelRight = right_cell_of_face[iface];
			if (cellLabelRight == -1) continue; //boundy face's right neighbor

			if ((cellLabelRight < 0)||(cellLabelRight >= nTotalCell)) {
				printf("Error: cellLabelRight = %d\n", cellLabelRight);
				exit(1);
			}
			int cell2faceOffset = offsetCell2Face[cellLabelLeft];
			cell2face[cell2faceOffset] = iface;
			//offsetCell2Face update for the next set of face
			offsetCell2Face[cellLabelLeft]++;
			cell2faceOffset = offsetCell2Face[cellLabelRight];
			cell2face[cell2faceOffset] = iface;
			offsetCell2Face[cellLabelRight]++;
		}		

		//reset offsetCell2Face
		printf("Before reset offsetCell2Face[%d]=%d, nTotalCellFaces = %d\n", nTotalCell-1, offsetCell2Face[nTotalCell-1], nTotalCellFaces);
		offsetCell2Face[0] = 0;
		for (int icell = 1; icell < nTotalCell; icell++){
			//offsetCell2Face[icell] = face_number_of_each_cell[icell - 1];
			offsetCell2Face[icell] = offsetCell2Face[icell-1] + face_number_of_each_cell[icell-1];
		}
		printf("After reset offsetCell2Face[%d]=%d, face_number_of_each_cell[%d] = %d, nTotalCellFaces = %d\n", nTotalCell-1, offsetCell2Face[nTotalCell-1], nTotalCell-1, face_number_of_each_cell[nTotalCell-1], nTotalCellFaces);
		//set cell2Cell
		//test begin
			//printf("right_cell_of_face[%d]=%d\n", 23221, right_cell_of_face[23221]);
		//test end
		for (int icell = 0; icell < nTotalCell; icell++){
			int sizeCellFace = face_number_of_each_cell[icell];
			int offsetCellFace = offsetCell2Face[icell];
			for (int iface = 0; iface < sizeCellFace; iface++){
				int faceID = cell2face[offsetCellFace + iface];
				int cellLabel = left_cell_of_face[faceID];
				if ((cellLabel < 0)||(cellLabel >= nTotalCell)) {
					printf("Error: cellLabel = %d\n", cellLabel);
					exit(1);
				}

				//test begin
					//if (icell == 0) printf("On the left, sizeCellFace =%d, iface = %d, faceID = %d, cellLabel = %d, rightCell = %d\n", sizeCellFace, iface, faceID, cellLabel, right_cell_of_face[faceID]);
				//test end

				if (cellLabel != icell) cell2cell[offsetCellFace + iface] = cellLabel;

				cellLabel = right_cell_of_face[faceID];
				//test begin
					//if (icell == 0) printf("On the right, faceID = %d, cellLabel = %d\n", faceID, cellLabel);
				//test end
				if (cellLabel == -1) continue; //boundy face's right neighbor
				if ((cellLabel < 0)||(cellLabel >= nTotalCell)) {
					printf("Error: cellLabel = %d\n", cellLabel);
					exit(1);
				}

				if (cellLabel != icell) cell2cell[offsetCellFace + iface] = cellLabel;
				
			}
		}

		//set bandwidth_of_each_cell
		computeBandwidth(nTotalCell, face_number_of_each_cell, offsetCell2Face, cell2cell);
		/*
		for (int icell = 0; icell < nTotalCell; icell++){
			int sizeCellFace = face_number_of_each_cell[icell];
			int offsetCellFace = offsetCell2Face[icell];

			//The max band between icell and neighbor cells
			int maxLocalBand = 0;
			for (int iface = 0; iface < sizeCellFace; iface++){
				int cellID = cell2cell[offsetCellFace+iface];
				int bandDiff = abs(icell - cellID);
				if (bandDiff > maxLocalBand) maxLocalBand = bandDiff;
				//test begin
					//if (icell == 0) printf("offsetCellFace = %d, iface = %d, cellID = %d, bandDiff = %d\n", offsetCellFace, iface, cellID, bandDiff);
				//test end
			}
			bandwidth_of_each_cell[icell] = maxLocalBand;
		}
		//statistics of bandwidth_of_each_cell
		int maxGlobalBand = 0; //maximum bandwidth
		int maxCellID = 0; //cellID of maximum bandwidth
		int bandDist[33]; //distribution of bandwidth from 1 to 32 and more
		for (int i = 0; i < 33; i++){
			bandDist[i] = 0;
		}

		for (int icell = 0; icell < nTotalCell; icell++){
			int bandCell = bandwidth_of_each_cell[icell];
			if (bandCell <= 0) {
				printf("Error: bandwidth<=0 on cell %d, bandwidth=%d\n", icell, bandCell);
				exit(1);
			}
			if (bandCell > maxGlobalBand) {
				maxGlobalBand = bandCell;
				maxCellID = icell;
			}
			if (bandCell > 32) bandDist[32]++;
			else bandDist[bandCell-1]++;
		}
		
		printf("maxGlobalBand = %d on cellID %d\n", maxGlobalBand, maxCellID);
		int nTotalDist = 0;
		for (int i = 0; i < 33; i++){
			nTotalDist += bandDist[i];
			printf("bandDist %d: %d in nTotalCells\n", i+1, bandDist[i]);
		}

		printf("nTotalDist = %d, nTotalCell = %d\n", nTotalDist, nTotalCell);
		*/
		delete []bandwidth_of_each_cell;		
	}
	
	void CallReOrderFaceLabel(Base_Grid_Conn * gconn){
		printf("Program is running in CallReOrderFaceLabel\n");
		cgsize_t nTotalFace = gconn->GetNTotalFace();
		cgsize_t nBoundFace = gconn->GetNBoundFace();
		cgsize_t nTotalCell = gconn->GetNTotalCell();
		printf("nTotalFace = %d, nBoundFace = %d, nTotalCell = %d\n", nTotalFace, nBoundFace, nTotalCell);

		vector<cgsize_t> &left_cell_of_face  = *gconn->GetLeftCellOfFace();
    		vector<cgsize_t> &right_cell_of_face = *gconn->GetRightCellOfFace();
		int * Oldleft_cell_of_face = new int[nTotalFace];
		int * Oldright_cell_of_face = new int[nTotalFace];
		for (int iFace = 0; iFace < nTotalFace; iFace++){
			Oldleft_cell_of_face[iFace] = left_cell_of_face[iFace];
			Oldright_cell_of_face[iFace] = right_cell_of_face[iFace];
			//add for test
				//if (right_cell_of_face[iFace] == -1) printf("right_cell_of_face on iFace %d is equal to -1\n");
			//add end
		}
		int * cell2face;
		int * leftRightCellFaces;//corresponding to cell2face, leftRightCellFaces record the cell is the left or right cell of one face. 

		//int * cell2cell;
		//offset of icell's cell faces in cell2face
		int * offsetCell2Face = new int[nTotalCell];
		int nTotalCellFaces = 0;

		int * face_number_of_each_cell = new int[nTotalCell];

		for (int icell = 0; icell < nTotalCell; icell++){
			face_number_of_each_cell[icell] = 0;
			offsetCell2Face[icell] = 0;
		}
		//set face_number_of_each_cell and nTotalCellFaces
		printf("set face_number_of_each_cell and nTotalCellFaces\n");
		SetFaceNumberOfEachCell(nTotalFace, nTotalCell, Oldleft_cell_of_face, Oldright_cell_of_face, face_number_of_each_cell, nTotalCellFaces);
		//set the size of cell2cell and cell2face
		printf("Before create of cell2cell nTotalCellFaces = %d\n", nTotalCellFaces);
		cell2face = new int[nTotalCellFaces];
		leftRightCellFaces = new int[nTotalCellFaces];
		//cell2cell = new int[nTotalCellFaces];

		//set offsetCell2Face
		SetOffSetCell2Face(nTotalCell, face_number_of_each_cell, offsetCell2Face);
		//set cell2face
		SetCell2Face(nTotalFace, nTotalCell, Oldleft_cell_of_face, Oldright_cell_of_face, offsetCell2Face, cell2face);
		//reset offsetCell2Face
		printf("Before reset offsetCell2Face[%d]=%d, nTotalCellFaces = %d\n", nTotalCell-1, offsetCell2Face[nTotalCell-1], nTotalCellFaces);
		//offsetCell2Face is modified due to update of cell2face
		//set offsetCell2Face again
		SetOffSetCell2Face(nTotalCell, face_number_of_each_cell, offsetCell2Face);
		printf("After reset offsetCell2Face[%d]=%d, face_number_of_each_cell[%d] = %d, nTotalCellFaces = %d\n", nTotalCell-1, offsetCell2Face[nTotalCell-1], nTotalCell-1, face_number_of_each_cell[nTotalCell-1], nTotalCellFaces);
		//check cell2face
		//rule 1: nBoundFace-1<cell2face<nTotalFace
		//rule 2: in cell2face, every value owns at least 1 item.
		printf("Checking cell2face ...\n");
		int * faceIDRecord = new int[nTotalFace-nBoundFace];
		for (int iFace = 0; iFace < (nTotalFace-nBoundFace); iFace++){
			faceIDRecord[iFace] = 0;
		}
		for (int cellID = 0; cellID < nTotalCell; cellID++){
			int numFaces = face_number_of_each_cell[cellID];
			int offset = offsetCell2Face[cellID];
			for (int iFace=0; iFace<numFaces; iFace++){
				int faceID = cell2face[offset+iFace];
				//rule 1
				if ((faceID < nBoundFace)||(faceID > nTotalFace)) {
					printf("Error: faceID %d on cellID %d and iFace %d is outside of nBoundFace %d and nTotalFace-1 %d\n", faceID, cellID, iFace, nBoundFace, nTotalFace-1);
					exit(1);
				}
				int remains = faceID - nBoundFace;
				faceIDRecord[remains]++;
			}
		}
		for (int iFace = 0; iFace < (nTotalFace-nBoundFace); iFace++){
			int numFaces = faceIDRecord[iFace];
			//rule 2
			if (numFaces < 1) {
				printf("Error: the number of faceID %d is %d\n", iFace+nBoundFace, numFaces);
				exit(1);
			}
		}
		
		//set leftRightCellFaces, which stores face
		//SetLeftRightCellFaces();
		int numBoundaryFaces = 0;
		for (int cellID = 0; cellID < nTotalCell; cellID++){
			int offset = offsetCell2Face[cellID];
			int numFaces = face_number_of_each_cell[cellID];
			for (int iFace = 0; iFace < numFaces; iFace++){
				int faceID = cell2face[offset+iFace];
				if (Oldleft_cell_of_face[faceID] == -1){
					printf("Error: on faceID, left_cell_of_face is -1\n", faceID);
					exit(1);
				}
				if (Oldright_cell_of_face[faceID] == -1) numBoundaryFaces++;
				//if (Oldright_cell_of_face[faceID] == -1) numBoundaryFaces++;
				if (Oldleft_cell_of_face[faceID] == cellID) leftRightCellFaces[offset+iFace] = 1;
				if (Oldright_cell_of_face[faceID] == cellID) leftRightCellFaces[offset+iFace] = -1;
				
			}
		}
		//if (numBoundaryFaces != nBoundFace) {
		//Because cell2face does not contain boudnary faces
		if (numBoundaryFaces != 0) {
			printf("Error: numBoundaryFaces=%d is not equal to 0\n", numBoundaryFaces);
			exit(1);
		}

		//Reorder faceID from nBoundFace to nTotalFace
		int * NewFaces = new int[nTotalFace]; //Old faceID and containerID
		//set the initial value of labelFace as nBoundFace-1
		int labelFace = nBoundFace-1;
		//Initialize NewFaces by -1, which means not changed.
		for (int faceID = 0; faceID < nTotalFace; faceID++){
			NewFaces[faceID] = -1;
		}
	     	#ifdef BSFCELLFACEORG
		//BSFCELLFACE means Broad Search First: faceID is distributed face by face in one cell with the ascending order. One cell starts before faceID distribution in one cell is finished. In fact, it is based on original faceID order. Before ReOrder, cell2face has been reordered with old faceID in ascending order.
			//Ascending cell2face for each cell
			printf("BSFCELLFACEORG: Ascending cell2face order...\n");
			for (int cellID = 0; cellID < nTotalCell; cellID++){
				int offset = offsetCell2Face[cellID];
				int numFaces = face_number_of_each_cell[cellID];
				int loopNum = numFaces - 1;
				//bubble sort
				for (int jFace = 0; jFace < numFaces-1; jFace++){
					for (int iFace = 0; iFace < loopNum; iFace++){
						int faceID = cell2face[offset+iFace];
						int faceIDNext = cell2face[offset+iFace+1];
						if (faceID > faceIDNext){
							cell2face[offset+iFace] = faceIDNext;
							cell2face[offset+iFace+1] = faceID;
							int temp = leftRightCellFaces[offset+iFace];
							leftRightCellFaces[offset+iFace] = leftRightCellFaces[offset+iFace+1];
							leftRightCellFaces[offset+iFace+1] = temp;
						}
					}

					loopNum--;

				}
			}
			printf("BSFCELLFACEORG: checking cell2face order...\n");
			for (int cellID = 0; cellID < nTotalCell; cellID++){
				int offset = offsetCell2Face[cellID];
				int numFaces = face_number_of_each_cell[cellID];
				for (int iFace = 0; iFace < numFaces-1; iFace++){
					int faceID = cell2face[offset+iFace];
					int faceIDNext = cell2face[offset+iFace+1];
					if (faceID > faceIDNext){
						printf("Error: faceID %d on cellID %d iFace %d is larger than faceIDNext %d\n", faceID, cellID, iFace, faceIDNext);
						exit(1);
					}
				}
			}
			printf("BSFCELLFACEORG: ReOrder cell labels ...\n");
			cellReorderByCell2FaceOrder(nTotalCell, nTotalFace, offsetCell2Face, face_number_of_each_cell, cell2face, nBoundFace, labelFace, NewFaces);
	     	#endif
	     	#ifdef BSFCELLFACELEFT
		//It is quite similar as BSFCELLFACEORG. Firstly, cell2face is reorded with old faceID with left_cell_of_face in ascending order. Then, set those faceID from nBoundFace to nTotalFace
			//Ascending cell2face for each cell
			printf("BSFCELLFACELEFT: Ascending cell2face order by left_cell_of_face ...\n");
			for (int cellID = 0; cellID < nTotalCell; cellID++){
				int offset = offsetCell2Face[cellID];
				int numFaces = face_number_of_each_cell[cellID];
				int loopNum = numFaces-1;
				//bubble sort
				for (int jFace = 0; jFace < numFaces-1; jFace++){
					for (int iFace = 0; iFace < loopNum; iFace++){
						int faceID = cell2face[offset+iFace];
						int leftCell = Oldleft_cell_of_face[faceID];
						int faceIDNext = cell2face[offset+iFace+1];
						int leftCellNext = Oldleft_cell_of_face[faceIDNext];
						if (leftCell > leftCellNext){
							cell2face[offset+iFace] = faceIDNext;
							cell2face[offset+iFace+1] = faceID;
							int temp = leftRightCellFaces[offset+iFace];
							leftRightCellFaces[offset+iFace] = leftRightCellFaces[offset+iFace+1];
							leftRightCellFaces[offset+iFace+1] = temp;
						}
					}
					loopNum--;
				}
			}
			printf("BSFCELLFACELEFT: checking cell2face order...\n");
			for (int cellID = 0; cellID < nTotalCell; cellID++){
				int offset = offsetCell2Face[cellID];
				int numFaces = face_number_of_each_cell[cellID];
				for (int iFace = 0; iFace < numFaces-1; iFace++){
					int faceID = cell2face[offset+iFace];
					int faceIDNext = cell2face[offset+iFace+1];
					int leftCell = Oldleft_cell_of_face[faceID];
					int leftCellNext = Oldleft_cell_of_face[faceIDNext];
					if (leftCell > leftCellNext){
						printf("Error: faceID %d on cellID %d iFace %d with leftCell %d is larger than leftCellNext %d\n", faceID, cellID, iFace, leftCell, leftCellNext);
						exit(1);
					}
				}
			}
			//Reorder cell label by faces in cell2face
			cellReorderByCell2FaceOrder(nTotalCell, nTotalFace, offsetCell2Face, face_number_of_each_cell, cell2face, nBoundFace, labelFace, NewFaces);
		#endif
	     	#ifdef BSFCELLFACERIGHT
		//It is quite similar as BSFCELLFACEORG. Firstly, cell2face is reorded with old faceID with right_cell_of_face in ascending order. Then, set those faceID from nBoundFace to nTotalFace
			//Ascending cell2face for each cell
			printf("BSFCELLFACERIGHT: Ascending cell2face order by right_cell_of_face ...\n");
			for (int cellID = 0; cellID < nTotalCell; cellID++){
				int offset = offsetCell2Face[cellID];
				int numFaces = face_number_of_each_cell[cellID];
				int loopNum = numFaces-1;
				//bubble sort
				for (int jFace = 0; jFace < numFaces-1; jFace++){
					for (int iFace = 0; iFace < loopNum; iFace++){
						int faceID = cell2face[offset+iFace];
						int rightCell = Oldright_cell_of_face[faceID];
						int faceIDNext = cell2face[offset+iFace+1];
						int rightCellNext = Oldright_cell_of_face[faceIDNext];
						if (rightCell > rightCellNext){
							cell2face[offset+iFace] = faceIDNext;
							cell2face[offset+iFace+1] = faceID;
							int temp = leftRightCellFaces[offset+iFace];
							leftRightCellFaces[offset+iFace] = leftRightCellFaces[offset+iFace+1];
							leftRightCellFaces[offset+iFace+1] = temp;
						}
					}
					loopNum--;
				}
			}
			printf("BSFCELLFACERIGHT: checking cell2face order...\n");
			for (int cellID = 0; cellID < nTotalCell; cellID++){
				int offset = offsetCell2Face[cellID];
				int numFaces = face_number_of_each_cell[cellID];
				for (int iFace = 0; iFace < numFaces-1; iFace++){
					int faceID = cell2face[offset+iFace];
					int faceIDNext = cell2face[offset+iFace+1];
					int rightCell = Oldright_cell_of_face[faceID];
					int rightCellNext = Oldright_cell_of_face[faceIDNext];
					if (rightCell > rightCellNext){
						printf("Error: faceID %d on cellID %d iFace %d with rightCell %d is larger than rightCellNext %d\n", faceID, cellID, iFace, rightCell, rightCellNext);
						exit(1);
					}
				}
			}
			//Reorder cell label by faces in cell2face
			cellReorderByCell2FaceOrder(nTotalCell, nTotalFace, offsetCell2Face, face_number_of_each_cell, cell2face, nBoundFace, labelFace, NewFaces);
		#endif
	     	#ifdef DEEPCELLFACE
		//DEEPCELLFACE means deep search first:in one cell, faceID is firstly distributed on a face with the smallest neighbor. Then, on the neighbour cell, faceID is distributed on a face with the smallest neighbor.
		//get the cell difference of left_cell_of_face and right_cell_of_face on faces
			int diffMax = 0;
			int diffMin = nTotalCell;
			//for (int faceID = 0; faceID < nTotalFace; faceID++){
			for (int faceID = nBoundFace; faceID < nTotalFace; faceID++){
				//boundary faces are not considered.
				//if ((Oldleft_cell_of_face[faceID] != - 1)&&(Oldright_cell_of_face[faceID] != -1)){
				int localDiff = Oldleft_cell_of_face[faceID] - Oldright_cell_of_face[faceID];
				if (diffMax < localDiff) diffMax = localDiff;
				if (diffMin > localDiff) diffMin = localDiff;
				//}
			}
			//In fact, it is also the maximumm bandwidth, which is a good way to compute the maximum bandwidth
			printf("diffMin = %d, diffMax = %d\n", diffMin, diffMax);
			int * diffNumMax = new int[diffMax];
			int * diffNumMin = new int[abs(diffMin)];
			//from 1 to diffMax
			for (int diffID = 0; diffID < diffMax; diffID++){
				diffNumMax[diffID] = 0;
			}
			//from -diffMin to -1
			for (int diffID = 0; diffID < abs(diffMin); diffID++){
				diffNumMin[diffID] = 0;
			}
			for (int faceID = nBoundFace; faceID < nTotalFace; faceID++){
				int localDiff = Oldleft_cell_of_face[faceID] - Oldright_cell_of_face[faceID];
				if (localDiff==0) printf("On faceID %d, localDiff=0\n", faceID);
				//if (localDiff >= 0) diffNumMax[localDiff]++;
				if (localDiff >= 0) diffNumMax[localDiff-1]++;
				//else diffNumMin[abs(localDiff)]++;
				else diffNumMin[abs(localDiff)-1]++;
			}
			//diff = 0 should not be repeated.
			//0~diffMax-1, -difMin~-1
			//0~diffMax-1, diffMax + abs(diffMin) -1 ~ diffMax
			int * offsetDiff = new int[diffMax + abs(diffMin) ];
			offsetDiff[0] = 0;
			//set diffMax firstly
			for (int diffID = 1; diffID < diffMax; diffID++){
				offsetDiff[diffID] = offsetDiff[diffID-1]+diffNumMax[diffID-1];
			}
			//set diffMin secondly
			//it is diffID = 0 for diffMin, which corresponds to -1
			offsetDiff[diffMax] = offsetDiff[diffMax-1]+diffNumMax[diffMax-1];
			for (int diffID = 1; diffID < abs(diffMin); diffID++){
				int offsetID = diffID + diffMax;
				//offsetDiff[offsetID] = offsetDiff[offsetID-1]+diffNumMin[diffID-1];
				offsetDiff[offsetID] = offsetDiff[offsetID-1]+diffNumMin[diffID-1];
			}
			for (int faceID = nBoundFace; faceID < nTotalFace; faceID++){
				int localDiff = Oldleft_cell_of_face[faceID] - Oldright_cell_of_face[faceID];
				if (localDiff >= 0) {
					//int newFaceID = offsetDiff[localDiff] + nBoundFace;
					//it must start from localDiff=1
					int newFaceID = offsetDiff[localDiff-1] + nBoundFace;
					NewFaces[newFaceID] = faceID;
					//offsetDiff[localDiff]++;
					offsetDiff[localDiff-1]++;
				}
				else {
					//int newFaceID = offsetDiff[abs(localDiff)+diffMax] + nBoundFace;
					//It may start from localDiff=-1
					int newFaceID = offsetDiff[abs(localDiff)+diffMax-1] + nBoundFace;
					NewFaces[newFaceID] = faceID;
					offsetDiff[abs(localDiff)+diffMax-1]++;
				}
			}

			//reset offsetDiff
			offsetDiff[0] = 0;
			//set diffMax firstly
			for (int diffID = 1; diffID < diffMax; diffID++){
				offsetDiff[diffID] = offsetDiff[diffID-1]+diffNumMax[diffID-1];
			}
			//set diffMin secondly
			//it is diffID = 0 for diffMin, which corresponds to -1
			offsetDiff[diffMax] = offsetDiff[diffMax-1]+diffNumMax[diffMax-1];
			for (int diffID = 1; diffID < abs(diffMin); diffID++){
				int offsetID = diffID + diffMax;
				//offsetDiff[offsetID] = offsetDiff[offsetID-1]+diffNumMin[diffID-1];
				offsetDiff[offsetID] = offsetDiff[offsetID-1]+diffNumMin[diffID-1];
			}

			//reorder NewFaces by left_cell_of_face ascending
			//diffID=0~diffMax corresponding to 1 to diffMax of left_cell_of_face - right_cell_of_face
			for (int diffID = 0; diffID < diffMax; diffID++){
				int numFaces = diffNumMax[diffID];
				int offset = offsetDiff[diffID];
				int startID = offset + nBoundFace;
				int endID = offset + nBoundFace + numFaces;
				int loopFaceID = endID - 1;
				//bubble sort newFaces
				for (int loopID = 0; loopID < numFaces-1; loopID++){
					for(int newFaceID = startID; newFaceID < loopFaceID; newFaceID++){
						int oldFaceID = NewFaces[newFaceID];
						int oldFaceIDNext = NewFaces[newFaceID+1];
						int leftCell = Oldleft_cell_of_face[oldFaceID];
						int leftCellNext = Oldleft_cell_of_face[oldFaceIDNext];
						if (leftCell > leftCellNext){
							NewFaces[newFaceID] = oldFaceIDNext;
							NewFaces[newFaceID+1] = oldFaceID;
						}
					}
					loopFaceID--;
				}
			}
			//check ascenidng order of NewFaces for diffMax
			printf("Checking diffMax ascending...\n");
			for (int diffID = 0; diffID < diffMax; diffID++){
				int numFaces = diffNumMax[diffID];
				int offset = offsetDiff[diffID];
				int startID = offset + nBoundFace;
				int endID = offset + nBoundFace + numFaces;
				for (int newFaceID = startID; newFaceID < endID-1; newFaceID++){
					int faceID = NewFaces[newFaceID];
					int faceIDNext = NewFaces[newFaceID+1];
					int leftCell = Oldleft_cell_of_face[faceID];
					int leftCellNext = Oldleft_cell_of_face[faceIDNext];
					if (leftCell > leftCellNext){
						printf("Error: on diffID %d, newFaceID = %d, faceID = %d, faceIDNext = %d, leftCell = %d, leftCellNext = %d\n", diffID, newFaceID, faceID, faceIDNext, leftCell, leftCellNext);
						exit(1);
					}

				}
			}
			//diffID=diffMax~diffMax+abs(dffMin) corresponding to -1 to diffMin of left_cell_of_face - right_cell_of_face
			for (int diffID = 0; diffID < abs(diffMin); diffID++){
				int diffMinID = diffMax + diffID;
				int numFaces = diffNumMin[diffID];
				int offset = offsetDiff[diffMinID];
				int startID = offset + nBoundFace;
				int endID = offset + nBoundFace + numFaces;
				int loopFaceID = endID - 1;
				//bubble sort newFaces
				for (int loopID = 0; loopID < numFaces-1; loopID++){
					for(int newFaceID = startID; newFaceID < loopFaceID; newFaceID++){
						int oldFaceID = NewFaces[newFaceID];
						int oldFaceIDNext = NewFaces[newFaceID+1];
						int leftCell = Oldleft_cell_of_face[oldFaceID];
						int leftCellNext = Oldleft_cell_of_face[oldFaceIDNext];
						if (leftCell > leftCellNext){
							NewFaces[newFaceID] = oldFaceIDNext;
							NewFaces[newFaceID+1] = oldFaceID;
						}
					}
					loopFaceID--;
				}
			}
			//check ascenidng order of NewFaces for diffMax
			printf("Checking diffMin ascending...\n");
			for (int diffID = 0; diffID < abs(diffMin); diffID++){
				int diffMinID = diffID + diffMax;
				int numFaces = diffNumMin[diffID];
				int offset = offsetDiff[diffMinID];
				int startID = offset + nBoundFace;
				int endID = offset + nBoundFace + numFaces;
				for (int newFaceID = startID; newFaceID < endID-1; newFaceID++){
					int faceID = NewFaces[newFaceID];
					int faceIDNext = NewFaces[newFaceID+1];
					int leftCell = Oldleft_cell_of_face[faceID];
					int leftCellNext = Oldleft_cell_of_face[faceIDNext];
					if (leftCell > leftCellNext){
						printf("Error: on diffID %d, newFaceID = %d, faceID = %d, faceIDNext = %d, leftCell = %d, leftCellNext = %d\n", diffID, newFaceID, faceID, faceIDNext, leftCell, leftCellNext);
						exit(1);
					}

				}
			}
		#endif
		//check reorder results
		//rule 1 Before nBoundFace, NewFaces=-1
		//rule 2 nBoundFace <= NewFaces < nTotalFace
		//rule 3 After nBoundFace, no elements are the same in NewFaces
		printf("Checking NewFaces which is ReOrder cell labels results...\n");
		for (int iFace = 0; iFace < nBoundFace; iFace++){
			//rule 1
			if (NewFaces[iFace]!=-1) {
				printf("Error in rule 1: on iFace %d, NewFaces is -1\n", iFace);
				exit(1);
			}
		}
		int * recordFaceID = new int[nTotalFace-nBoundFace];
		for (int iFace = 0; iFace < nTotalFace-nBoundFace; iFace++){
			recordFaceID[iFace] = 0;
		}
		for (int iFace = nBoundFace; iFace < nTotalFace-nBoundFace; iFace++){
			int faceID = NewFaces[iFace];
			//rule 2
			if ((faceID < nBoundFace)||(faceID>(nTotalFace-1))){
				printf("Error in rule 2: on iFace %d of NewFaces, faceID = %d is outside of nBoundFace %d and nTotalFace - 1 %d\n", iFace, faceID, nBoundFace, nTotalFace-1);
				exit(1);
			}
			int remains = (nTotalFace-1) - faceID;
			recordFaceID[remains]++;
		}
		//rule 3
		for (int iFace = 0; iFace < nTotalFace-nBoundFace; iFace++){
			if (recordFaceID[iFace] > 1) {
				printf("Error: more than two faces own the same faceID %d\n", nTotalFace - 1 - iFace);
				exit(1);
			}
		}
		
		//set NewFaces bwtween 0 and nBoundFace-1
		for (int iFace = 0; iFace < nBoundFace; iFace++){
			NewFaces[iFace] = iFace;
		}
	
		//create reflectFaces for recording relationship of containerID and old cellID
		int * reflectFaces = new int[nTotalFace];
		for (int oldFaceID = 0; oldFaceID < nTotalFace; oldFaceID++){
			int containerID = NewFaces[oldFaceID];
			reflectFaces[containerID] = oldFaceID;
		}
		//check reflectFaceID
		//rule 1: 0~nBoundFace, reflectFaces[iFace]=iFace
		//rule 2: nBoundFace~nTotalFace, nBoundFace<=reflectFaces<nTotalFace
		//rule 3: 0~nTotalFace, reflectFaces does not own the same old faceID
		printf("Checking reflectFaceID  results...\n");
		for (int iFace = 0; iFace < nBoundFace; iFace++){
			//rule 1
			if (reflectFaces[iFace] != iFace) {
				printf("Error: on containerID %d, is not the same faceID\n", iFace);
				exit(1);
			}
		}
		for (int iFace = 0; iFace < nTotalFace-nBoundFace; iFace++){
			recordFaceID[iFace] = 0;
		}
		for (int iFace = nBoundFace; iFace < nTotalFace; iFace++){
			int oldFaceID = reflectFaces[iFace];
			if ((oldFaceID < nBoundFace)||(oldFaceID > nTotalFace)) {
				//rule 2
				printf("Error: on containerID %d, oldFaceID %d is outside nBoundFace and nTotalFace\n", iFace, oldFaceID);
				exit(1);
			}
			int remains = oldFaceID - nBoundFace;
			recordFaceID[remains]++;
		}
		for (int iFace = 0; iFace < nTotalFace-nBoundFace; iFace++){
			if (recordFaceID[iFace] != 1) {
				//rule 3
				printf("Error: on iFace %d, face number %d is not 1\n", iFace+nBoundFace, recordFaceID[iFace]);
				exit(1);
			}
		}
		
		//swap method for updating node_number_of_each_face, face2node,
		//left_cell_of_face and right_cell_of_face
		vector< vector<cgsize_t> > & face2node = *gconn->GetFace2Node();
		vector<int> & node_number_of_each_face = *gconn->GetNodeNumberOfEachFace();
				
		int * OrgFace = new int[nTotalFace]; //containerID and old faceID reflection
		int * SwapFace = new int[nTotalFace];//old faceID and containerID reflection
		for (int iFace = 0; iFace < nTotalFace; iFace++){
			OrgFace[iFace] = iFace;
			SwapFace[iFace] = iFace;
		}

		for (int faceID = 0; faceID < nTotalFace; faceID++){
			if (OrgFace[faceID] == reflectFaces[faceID]) continue;
			int targetID = reflectFaces[faceID];
			int currentID = OrgFace[faceID];
			int swapID = SwapFace[targetID];
			//update node_number_of_each_face on faceID and swapID
			int tempNumNode = node_number_of_each_face[swapID];
			node_number_of_each_face[swapID] = node_number_of_each_face[faceID];
			node_number_of_each_face[faceID] = tempNumNode;
			//update face2node on faceID and swapID
			vector<cgsize_t> tempFace = face2node[swapID];
			face2node[swapID] = face2node[faceID];
			face2node[faceID] = tempFace;
			//update left_cell_of_face on faceID and swapID
			int tempCell = left_cell_of_face[swapID];
			left_cell_of_face[swapID] = left_cell_of_face[faceID];
			left_cell_of_face[faceID] = tempCell;
			//update right_cell_of_face on faceID and swapID
			tempCell = right_cell_of_face[swapID];
			right_cell_of_face[swapID] = right_cell_of_face[faceID];
			right_cell_of_face[faceID] = tempCell;
			//update OrgFace on faceID and swapID
			OrgFace[faceID] = targetID;
			OrgFace[swapID] = currentID;
			//update swapID on targetID and currentID
			SwapFace[targetID] = faceID;
			SwapFace[currentID] = swapID;
		}
		//check OrgFace and SwapFace
		printf("Check OrgFace and SwapFace ...\n");
		for (int iFace = 0; iFace < nTotalFace; iFace++){
			if (OrgFace[iFace] != reflectFaces[iFace]){
				printf("Error: on iFace %d, OrgFace %d is not equal to reflectFaces %d\n", iFace, OrgFace[iFace], reflectFaces[iFace]);
				exit(1);
			}

			if (SwapFace[iFace] != NewFaces[iFace]){
				printf("Error:  on iFace %d, SwapFace %d is not equal to NewFaces %d\n", iFace, SwapFace[iFace], NewFaces[iFace]);
				exit(1);
			}
		}
		
	}

	void CallReOrderCellLabel(Base_Grid_Conn * gconn){
		printf("Program is running in CallReOrderCellLabel\n");
		cgsize_t nTotalFace = gconn->GetNTotalFace();
		cgsize_t nBoundFace = gconn->GetNBoundFace();
		cgsize_t nTotalCell = gconn->GetNTotalCell();
		printf("nTotalFace = %d, nBoundFace = %d, nTotalCell = %d\n", nTotalFace, nBoundFace, nTotalCell);

		vector<cgsize_t> &left_cell_of_face  = *gconn->GetLeftCellOfFace();
    		vector<cgsize_t> &right_cell_of_face = *gconn->GetRightCellOfFace();
		int * Oldleft_cell_of_face = new int[nTotalFace];
		int * Oldright_cell_of_face = new int[nTotalFace];
		for (int iFace = 0; iFace < nTotalFace; iFace++){
			Oldleft_cell_of_face[iFace] = left_cell_of_face[iFace];
			Oldright_cell_of_face[iFace] = right_cell_of_face[iFace];
		}
		/*
		//Just a test
			printf("iface = %d, leftCell = %d, rightCell = %d\n", 451581, left_cell_of_face[451581], right_cell_of_face[451581]);
			printf("iface = %d, leftCell = %d, rightCell = %d\n", 451582, left_cell_of_face[451582], right_cell_of_face[451582]);
		//test end
		*/
		int * cell2face;

		int * cell2cell;
		//offset of icell's cell faces in cell2face
		int * offsetCell2Face = new int[nTotalCell];
		int nTotalCellFaces = 0;

		int * face_number_of_each_cell = new int[nTotalCell];

		for (int icell = 0; icell < nTotalCell; icell++){
			face_number_of_each_cell[icell] = 0;
			offsetCell2Face[icell] = 0;
		}

		//set face_number_of_each_cell and nTotalCellFaces
		printf("set face_number_of_each_cell and nTotalCellFaces\n");
		SetFaceNumberOfEachCell(nTotalFace, nTotalCell, Oldleft_cell_of_face, Oldright_cell_of_face, face_number_of_each_cell, nTotalCellFaces);
		/*
		for (int iface = 0; iface < nTotalFace; iface++){
			//for left_cell_of_face
			int cellLabelLeft = left_cell_of_face[iface];
			if ((cellLabelLeft < 0)||(cellLabelLeft >= nTotalCell)) {
				printf("Error: cellLabelLeft = %d\n", cellLabelLeft);
				exit(1);
			}
			int cellLabelRight = right_cell_of_face[iface];
			//Boundary faces are not considered. In fact, on boundary face, cellLabelLeft is cell itself. cellLabelRight is -1, which is not necessary stored in cell2cell.
			if (cellLabelRight == -1) continue; //boundy face's right neighbor

			if ((cellLabelRight < 0)||(cellLabelRight >= nTotalCell)) {
				printf("Error: cellLabelRight = %d\n", cellLabelRight);
				exit(1);
			}
			face_number_of_each_cell[cellLabelLeft]++;
			nTotalCellFaces++;
			face_number_of_each_cell[cellLabelRight]++;
			nTotalCellFaces++;
		}		
		*/
		//set the size of cell2cell and cell2face
		printf("Before create of cell2cell nTotalCellFaces = %d\n", nTotalCellFaces);
		cell2face = new int[nTotalCellFaces];
		cell2cell = new int[nTotalCellFaces];

		//set offsetCell2Face
		SetOffSetCell2Face(nTotalCell, face_number_of_each_cell, offsetCell2Face);
		/*
		offsetCell2Face[0] = 0;
		for (int icell = 1; icell < nTotalCell; icell++){
			offsetCell2Face[icell] = offsetCell2Face[icell-1] + face_number_of_each_cell[icell-1];
		}
		*/
		//set cell2face
		SetCell2Face(nTotalFace, nTotalCell, Oldleft_cell_of_face, Oldright_cell_of_face, offsetCell2Face, cell2face);
		/*
		for (int iface = 0; iface < nTotalFace; iface++){
			//for left_cell_of_face
			int cellLabelLeft = left_cell_of_face[iface];
			if ((cellLabelLeft < 0)||(cellLabelLeft >= nTotalCell)) {
				printf("Error: cellLabelLeft = %d\n", cellLabelLeft);
				exit(1);
			}

			int cellLabelRight = right_cell_of_face[iface];
			//boundary face is not considered.
			if (cellLabelRight == -1) continue; //boundy face's right neighbor

			if ((cellLabelRight < 0)||(cellLabelRight >= nTotalCell)) {
				printf("Error: cellLabelRight = %d\n", cellLabelRight);
				exit(1);
			}
			int cell2faceOffset = offsetCell2Face[cellLabelLeft];
			cell2face[cell2faceOffset] = iface;
			//offsetCell2Face update for the next set of face
			offsetCell2Face[cellLabelLeft]++;
			cell2faceOffset = offsetCell2Face[cellLabelRight];
			cell2face[cell2faceOffset] = iface;
			offsetCell2Face[cellLabelRight]++;
		}		
		*/
		//reset offsetCell2Face
		printf("Before reset offsetCell2Face[%d]=%d, nTotalCellFaces = %d\n", nTotalCell-1, offsetCell2Face[nTotalCell-1], nTotalCellFaces);
		//offsetCell2Face is modified due to update of cell2face
		//set offsetCell2Face again
		SetOffSetCell2Face(nTotalCell, face_number_of_each_cell, offsetCell2Face);
		/*
		offsetCell2Face[0] = 0;
		for (int icell = 1; icell < nTotalCell; icell++){
			//offsetCell2Face[icell] = face_number_of_each_cell[icell - 1];
			offsetCell2Face[icell] = offsetCell2Face[icell-1] + face_number_of_each_cell[icell-1];
		}
		*/
		printf("After reset offsetCell2Face[%d]=%d, face_number_of_each_cell[%d] = %d, nTotalCellFaces = %d\n", nTotalCell-1, offsetCell2Face[nTotalCell-1], nTotalCell-1, face_number_of_each_cell[nTotalCell-1], nTotalCellFaces);

		//set cell2Cell
		SetCell2Cell(nTotalCell, face_number_of_each_cell, offsetCell2Face, cell2face, Oldleft_cell_of_face, Oldright_cell_of_face, cell2cell);
		/*
		for (int icell = 0; icell < nTotalCell; icell++){
			int sizeCellFace = face_number_of_each_cell[icell];
			int offsetCellFace = offsetCell2Face[icell];
			//for test
				//if (257976 == icell) printf("sizeCellFace = %d in icell %d\n", sizeCellFace, icell);
			//test end
			for (int iface = 0; iface < sizeCellFace; iface++){
				int faceID = cell2face[offsetCellFace + iface];
				//just a test
				//if (257976 == icell) printf("icell = %d, faceID = %d, leftCell = %d, rightCell = %d\n", icell, faceID, left_cell_of_face[faceID], right_cell_of_face[faceID]);
				//test end
				int cellLabel = left_cell_of_face[faceID];
				if ((cellLabel < 0)||(cellLabel >= nTotalCell)) {
					printf("Error: cellLabel = %d\n", cellLabel);
					exit(1);
				}

				if (cellLabel != icell) cell2cell[offsetCellFace + iface] = cellLabel;

				cellLabel = right_cell_of_face[faceID];

				if (cellLabel == -1) continue; //boundy face's right neighbor
				if ((cellLabel < 0)||(cellLabel >= nTotalCell)) {
					printf("Error: cellLabel = %d\n", cellLabel);
					exit(1);
				}

				if (cellLabel != icell) cell2cell[offsetCellFace + iface] = cellLabel;
				
			}
		}
		*/
		//check cell2Cell, to see some rules:
		//rule 1 one cell's neighbor cells should not contain itself.
		//rule 2 0<= one cell's neighbor cells < nTotalCell
		//rule 3 one cell's neighbor cells should be different.
		//It should be noted that for O type grid, break 3 may be broken.
		CheckCell2Cell(nTotalCell, face_number_of_each_cell, offsetCell2Face, cell2cell);
		//By above process, cell2cell is established, which stores one cell's neighbour cells' labels. offsetCell2Face can be used by access of cell2cell. face_number_of_each_cell stores face number owned by one cell, which can also be used for neighbour cell number owned by one cell. Boundary faces are not considered in both cell2cell, offsetCell2Face and face_number_of_each_cell.
		//For reording cell index, cell2cell with offsetCell2Face can supply cell connection information. face_number_of_each_cell can supply degree of cell. 
		int maxDegree = 0;
		for (int iCell = 0; iCell < nTotalCell; iCell++){
			int degreeCell = face_number_of_each_cell[iCell];
			if (maxDegree < degreeCell) maxDegree = degreeCell;
		}
		printf("The maximum degree of cell is %d\n", maxDegree);
		//establish containers
		int * ReOrderCell = new int[nTotalCell]; //reorder cell labels
		int * NewCell = new int[nTotalCell]; //reorder cell labels
		//int maxSizeLevel = maxDegree * maxDegree; //maximum number of one level
		//It may be very large, not just maxDegree*maxDegree
		//int maxSizeLevel = 5000; //maximum number of one level
		//The worst condition
		int maxSizeLevel = nTotalCell; //maximum number of one level
		int * cellIDSortCur = new int[maxSizeLevel];
		int numSortCur = 0;
		int * cellIDSortNext = new int[maxSizeLevel];
		int numSortNext = 0;
		//Initialization of ReOrderCell by -1
		for (int iCell = 0; iCell < nTotalCell; iCell++) ReOrderCell[iCell] = -1;
		//Initialization of NewCell by -1
		for (int iCell = 0; iCell < nTotalCell; iCell++) NewCell[iCell] = -1;
		//Initialization of cellIDSortCur and cellIDSortNext by -1
		for (int iCell = 0; iCell < maxSizeLevel; iCell++) {
			cellIDSortCur[iCell] = -1;
			cellIDSortNext[iCell] = -1;
		}
		//Just a test
		//printf("Test cellID %d Before ReOrdering\n", 257976);
			//ReOrderCell[257976] = 0;	
			//GetCurrentLevelForReOrder(257976, face_number_of_each_cell, cell2cell, offsetCell2Face, ReOrderCell, numSortCur, cellIDSortCur);
			//TestCellIDSortCurDegree(numSortCur, cellIDSortCur, face_number_of_each_cell);
		//Test end
		//Searching the root vertex
		printf("Searching the root vertex ...\n");
		//start from vetex 0	
		int rStart = 0;
		printf("At the beginning, rStart = %d\n", rStart);
		//Here, ReOrderCell is used for marking a cell exsiting in root level structure
		ReOrderCell[rStart] = 1;
		int lCur = 0;
		int lOld = -1;
		//Get current level of vertex r
		GetCurrentLevel(rStart, face_number_of_each_cell, cell2cell, offsetCell2Face, ReOrderCell, numSortCur, cellIDSortCur);	
		/*
		//Test for current level of rStart = 0
		printf("numSortCur = %d\n", numSortCur);
		for (int iSort=0; iSort<numSortCur; iSort++){
			printf("iSort = %d, cellID = %d, cell2cell=%d\n", iSort, cellIDSortCur[iSort], cell2cell[iSort]);
		}
		*/
		//Get the next level of vertex r
		GetNextLevel(numSortCur, cellIDSortCur, face_number_of_each_cell, cell2cell, offsetCell2Face, ReOrderCell, numSortNext, cellIDSortNext);
		/*
		//Test for next level of rStart=0
		printf("numSortNext=%d\n", numSortNext);
		for (int iSort=0; iSort < numSortNext; iSort++){
			printf("iSort = %d of numSortNext %d, cellID = %d\n", iSort, numSortNext, cellIDSortNext[iSort]);
		}
		*/
		while (lCur > lOld){
			lOld = lCur;
			lCur = 0;
			while (numSortNext > 0){
				lCur++;
				/*
				//Test numSortCur and cellIDSortCur before SetCellIDSortCurByCellIDSortNext
				printf("Test numSortCur and cellIDSortCur before SetCellIDSortCurByCellIDSortNext\n");
				TestCellSort(numSortCur, cellIDSortCur);
				*/
				
				ResetCellSort(numSortCur, cellIDSortCur);
				SetCellIDSortCurByCellIDSortNext(numSortNext, cellIDSortNext, numSortCur, cellIDSortCur);
				ResetCellSort(numSortNext, cellIDSortNext);
				/*
				//Test numSortCur and cellIDSortCur after SetCellIDSortCurByCellIDSortNext
				printf("Test numSortCur and cellIDSortCur after SetCellIDSortCurByCellIDSortNext\n");
				TestCellSort(numSortCur, cellIDSortCur);
				*/
				GetNextLevel(numSortCur, cellIDSortCur, face_number_of_each_cell, cell2cell, offsetCell2Face, ReOrderCell, numSortNext, cellIDSortNext);

			}


			if (lCur > lOld){
				/*
				//Test numSortCur and cellIDSortCur after while (numSortNext > 0)
				printf("Test numSortCur and cellIDSortCur after while (numSortNext > 0)\n");
				TestCellSort(numSortCur, cellIDSortCur);
				*/
				//get the minimum degree in cellIDSortCur, start the next iteration
				rStart = FindMinCellIDSortCur(numSortCur, cellIDSortCur, face_number_of_each_cell);
				//printf("rStart = %d in numSortCur %d\n", rStart, numSortCur);
				//reset R by -1
				for (int iCell = 0; iCell < nTotalCell; iCell++) ReOrderCell[iCell] = -1;
				//reset numSortCur by -1
				ResetCellSort(numSortCur, cellIDSortCur);
				//reset numSortNext by -1
				ResetCellSort(numSortNext, cellIDSortNext);

				ReOrderCell[rStart] = 1;

				GetCurrentLevel(rStart, face_number_of_each_cell, cell2cell, offsetCell2Face, ReOrderCell, numSortCur, cellIDSortCur);	

				GetNextLevel(numSortCur, cellIDSortCur, face_number_of_each_cell, cell2cell, offsetCell2Face, ReOrderCell, numSortNext, cellIDSortNext);
					
				
			}
			printf("lCur = %d, lOld = %d, numSortNext=%d, rStart = %d\n", lCur, lOld, numSortNext ,rStart);
		}

		printf("Reorder cell index ...\n");
		printf("Start from %d\n", rStart);
		//reset ReOrderCell by -1
		for (int iCell = 0; iCell < nTotalCell; iCell++) ReOrderCell[iCell] = -1;
		//reset numSortCur by -1
		ResetCellSort(numSortCur, cellIDSortCur);

		//reset numSortNext by -1
		ResetCellSort(numSortNext, cellIDSortNext);
	
		//set rStart as 0, and numVertex as 1
		ReOrderCell[rStart] = 1;	
		int numVertex = 0;
		NewCell[rStart] = numVertex;	
		int nLevel = 0;
		
		//Get current level of rStart
		//GetCurrentLevelForReOrder(rStart, face_number_of_each_cell, cell2cell, offsetCell2Face, ReOrderCell, numSortCur, cellIDSortCur);
		GetCurrentLevel(rStart, face_number_of_each_cell, cell2cell, offsetCell2Face, ReOrderCell, numSortCur, cellIDSortCur);

		while(numSortCur > 0){
			printf("nLevel = %d, numSortCur = %d\n", nLevel, numSortCur);
			//printf("Before AscendingCellIDSortCurByDegree\n");
			//TestCellIDSortCurDegree(numSortCur, cellIDSortCur, face_number_of_each_cell);
			AscendingCellIDSortCurByDegree(numSortCur, face_number_of_each_cell, cellIDSortCur);	
			//printf("After AscendingCellIDSortCurByDegree\n");
			//TestCellIDSortCurDegree(numSortCur, cellIDSortCur, face_number_of_each_cell);
			//printf("Before ReOrderLevelCur, numVertex = %d\n", numVertex);
			//ReOrderLevelCur(numSortCur, cellIDSortCur, numVertex, ReOrderCell);
			ReOrderLevelCur(numSortCur, cellIDSortCur, numVertex, NewCell);
			//printf("After ReOrderLevelCur, numVertex = %d\n", numVertex);

			//GetNextLevelForReOrder(numSortCur, cellIDSortCur, face_number_of_each_cell, cell2cell, offsetCell2Face, ReOrderCell, numSortNext, cellIDSortNext);
			GetNextLevel(numSortCur, cellIDSortCur, face_number_of_each_cell, cell2cell, offsetCell2Face, ReOrderCell, numSortNext, cellIDSortNext);

			ResetCellSort(numSortCur, cellIDSortCur);
			SetCellIDSortCurByCellIDSortNext(numSortNext, cellIDSortNext, numSortCur, cellIDSortCur);
			ResetCellSort(numSortNext, cellIDSortNext);
			nLevel++;

		}
		printf("Finish cell reorder by RCM scheme...\n");
		printf("nToatalCell = %d, numVertex = %d\n", nTotalCell, numVertex);
		//reverse NewCell
		//Not sure
		for (int iCell = 0; iCell < nTotalCell; iCell++){
			NewCell[iCell] = nTotalCell - 1 - NewCell[iCell];
		}
		//check rules for cell index
		CheckNewCell(nTotalCell, NewCell);
		//reflectCell is a reverse reflect from NewCell to 
		int * reflectCell = new int[nTotalCell];
		for (int iCell = 0; iCell < nTotalCell; iCell++){
			int cellID = NewCell[iCell];
			reflectCell[cellID] = iCell;
		}
		//define new face_number_of_each_cell and new offsetCell2Face
		int * NewFace_number_of_each_cell = new int[nTotalCell];
		int * NewOffsetCell2Face = new int[nTotalCell];
		//set NewFace_number_of_each_cell
		for (int cellID = 0; cellID < nTotalCell; cellID++){
			int reflectCellID = reflectCell[cellID];
			NewFace_number_of_each_cell[cellID] = face_number_of_each_cell[reflectCellID];
		}

		SetOffSetCell2Face(nTotalCell, NewFace_number_of_each_cell, NewOffsetCell2Face);
		//define new left_cell_of_face and right_cell_of_face
		//define new cell2cell and cell2face
		printf("Before create of NewCell2Cell, nTotalCellFaces = %d\n", nTotalCellFaces);
		//int * NewCell2Face = new int[nTotalCellFaces];
		int * NewCell2Cell = new int[nTotalCellFaces];
		//set NewCell2Cell
		for (int cellID = 0; cellID < nTotalCell; cellID++){
			int reflectCellID = reflectCell[cellID];
			int numCellCellsOld = face_number_of_each_cell[reflectCellID];
			int numCellCellsNew = NewFace_number_of_each_cell[cellID];
			if (numCellCellsOld != numCellCellsNew) {
				printf("Error: on cellID %d, numCellCellsNew %d is not equal to numCellCellsOld %d on reflectCellID %d\n", cellID, numCellCellsNew, numCellCellsOld, reflectCellID);
				exit(1);
			}
			int offsetNew = NewOffsetCell2Face[cellID];
			int offsetOld = offsetCell2Face[reflectCellID];
			for (int iCell = 0; iCell < numCellCellsNew; iCell++){
				int newCellID = NewCell[cell2cell[offsetOld+iCell]];
				//NewCell2Cell[offsetNew+iCell] = cell2cell[offsetOld+iCell];
				NewCell2Cell[offsetNew+iCell] = newCellID;
			}
		}
		//set new cell2Cell
		//SetCell2Cell(nTotalCell, Newface_number_of_each_cell, NewoffsetCell2Face, NewCell2Face, Newleft_cell_of_face, Newright_cell_of_face, NewCell2Cell);
	
		//check NewCell2Cell
		CheckCell2Cell(nTotalCell, NewFace_number_of_each_cell, NewOffsetCell2Face, NewCell2Cell);
		//calculate bandwidth again
			
		computeBandwidth(nTotalCell, NewFace_number_of_each_cell, NewOffsetCell2Face, NewCell2Cell);
		//According to cell reorder, Update face_number_of_each_face by swap
		//Input: face_nuber_of_each_cell and reflectCell
		SwapFaceNumberOfEachFace(nTotalCell, face_number_of_each_cell, reflectCell, NewCell, NewFace_number_of_each_cell);
		//Update left_cell_of_face and right_cell_of_face
		for (int iFace = 0; iFace < nTotalFace; iFace++){
			//for left_cell_of_face
			int OldCellID = left_cell_of_face[iFace];
			//not bounday faces, treatment is requierd
			if (OldCellID != -1){
				int NewCellID = NewCell[OldCellID];
				left_cell_of_face[iFace] = NewCellID;
			}
			//for right_cell_of_face
			OldCellID = right_cell_of_face[iFace];
			if (OldCellID != -1){
				int NewCellID = NewCell[OldCellID];
				right_cell_of_face[iFace] = NewCellID;

			}
		}
		//update cell2node and node_number_of_each_cell by swap method
		vector< vector<cgsize_t> > &cell2node = *gconn->GetCell2Node();
		vector<int> & node_number_of_each_cell = *gconn->GetNodeNumberOfEachCell();
		int * OrgCell = new int[nTotalCell];//Contianer index and old cellID reflection
		int * SwapCell = new int[nTotalCell]; //Old cellID and container index reflection
		for (int cellID = 0; cellID < nTotalCell; cellID++){
			//Before swap, OrgCell is equal to SwapCell
			OrgCell[cellID] = cellID; 
			SwapCell[cellID] = cellID;
		}

		for (int cellID = 0; cellID < nTotalCell; cellID++){
			if (OrgCell[cellID] == reflectCell[cellID]) continue;
			int targetID = reflectCell[cellID];
			int currentID = OrgCell[cellID];
			int SwapID = SwapCell[targetID];
			//update correspondings of cellID and SwapID for node_number_of_each_cell
			int tempNodeNumber = node_number_of_each_cell[SwapID];
			node_number_of_each_cell[SwapID] = node_number_of_each_cell[cellID];
			node_number_of_each_cell[cellID] = tempNodeNumber;
			//update correspondings of cellID and SwapID for cell2node
			vector<cgsize_t> tempCell = cell2node[SwapID];
			cell2node[SwapID] = cell2node[cellID];
			cell2node[cellID] = tempCell;
			//update correspondings of cellID and SwapID for OrgCell
			OrgCell[cellID] = targetID;
			OrgCell[SwapID] = currentID;
			//update correspondings of targerID and currentID for SwapID
			SwapCell[targetID] = cellID;
			SwapCell[currentID] = SwapID;
		}
		

	}

	void GetCurrentLevel(int rStart, const int * cell_number_of_each_cell, const int * cell2cell, const int * offsetCell2Cell, int * ReOrderCell, int & numSortCur, int * cellIDSortCur){
		numSortCur = 0;
		int numCellCells = cell_number_of_each_cell[rStart];
		int offset = offsetCell2Cell[rStart];
		for (int iCell = 0; iCell < numCellCells; iCell++){
			int cellID = cell2cell[offset + iCell];
			//cellID is not in the root level structure
			if (ReOrderCell[cellID] == -1){
				cellIDSortCur[numSortCur] = cellID;
				numSortCur++;
				ReOrderCell[cellID] = 1;
			}
		}
	}

	void GetNextLevel(const int & numSortCur, const int * cellIDSortCur, const int * cell_number_of_each_cell, const int * cell2cell, const int * offsetCell2Cell, int * ReOrderCell, int & numSortNext, int * cellIDSortNext){
		numSortNext = 0;
		for (int iCell = 0; iCell < numSortCur; iCell++){
			int cellCellID = cellIDSortCur[iCell];
			int numCellCells = cell_number_of_each_cell[cellCellID];
			int offset = offsetCell2Cell[cellCellID];
			for (int iCellCell = 0; iCellCell < numCellCells; iCellCell++){
				int cellID = cell2cell[offset + iCellCell];
				if (ReOrderCell[cellID] == -1){
					cellIDSortNext[numSortNext] = cellID;
					numSortNext++;
					ReOrderCell[cellID] = 1;
				}
			}
		}
	}

	void SetCellIDSortCurByCellIDSortNext(const int & numSortNext, const int * cellIDSortNext, int & numSortCur, int * cellIDSortCur){
		numSortCur = numSortNext;
		for (int iSort=0; iSort < numSortNext; iSort++){
			cellIDSortCur[iSort] = cellIDSortNext[iSort];
		}
	}

	void ResetCellSort(int & numSort, int * cellIDSort){
		for (int iSort=0; iSort<numSort; iSort++){
			cellIDSort[iSort] = -1;
		}
		numSort = 0;
	}

	void TestCellSort(const int numSort, const int * cellIDSort){
		for (int iSort = 0; iSort<numSort; iSort++){
			printf("iSort = %d of numSort %d, cellID = %d\n", iSort, numSort, cellIDSort[iSort]);
		}
	}

	int FindMinCellIDSortCur(const int numSortCur, const int * cellIDSortCur, const int * cell_number_of_each_cell){
		int minCellID = cellIDSortCur[0];
		int minDegree = cell_number_of_each_cell[minCellID];
		for (int iSort = 1; iSort < numSortCur; iSort++){
			int cellID = cellIDSortCur[iSort];
			if (cell_number_of_each_cell[iSort] < minDegree){
				minDegree = cell_number_of_each_cell[iSort];
				minCellID = cellID;
			}

		}

		return minCellID;
	}

	void GetCurrentLevelForReOrder(int rStart, const int * cell_number_of_each_cell, const int * cell2cell, const int * offsetCell2Cell, int * ReOrderCell, int & numSortCur, int * cellIDSortCur){
		numSortCur = 0;
		int numCellCells = cell_number_of_each_cell[rStart];
		//printf("numCellCells = %d of cellID = %d\n", numCellCells, rStart);
		int offset = offsetCell2Cell[rStart];
		for (int iCell = 0; iCell < numCellCells; iCell++){
			int cellID = cell2cell[offset + iCell];
			//printf("iCell = %d, cellID = %d\n", iCell, cellID);
			//cellID is not in the root level structure
			if (ReOrderCell[cellID] == -1){
				cellIDSortCur[numSortCur] = cellID;
				numSortCur++;
			}
		}
	}
	
	void AscendingCellIDSortCurByDegree(const int numSortCur, const int * cell_number_of_each_cell, int * cellIDSortCur){
		int numVertexInSort = numSortCur;
		//At least 2 vertices for sort
		while(numVertexInSort > 1){
			for (int iSort=0; iSort < numVertexInSort-1; iSort++){
				int cellID = cellIDSortCur[iSort];
				int cellIDNext = cellIDSortCur[iSort+1];
				int degree = cell_number_of_each_cell[cellID];
				int degreeNext = cell_number_of_each_cell[cellIDNext];
				if (degree > degreeNext){
					//exchange vertex order
					cellIDSortCur[iSort] = cellIDNext;		
					cellIDSortCur[iSort+1] = cellID;
				}
			}
			numVertexInSort--;
		}
	}

	void TestCellIDSortCurDegree(const int numSortCur, const int * cellIDSortCur, const int * cell_number_of_each_cell){
		//Test Ascending order of vertex 0
		for (int iSort = 0; iSort < numSortCur; iSort++){
			printf("iSort = %d, cellID = %d, degree = %d\n", iSort, cellIDSortCur[iSort], cell_number_of_each_cell[cellIDSortCur[iSort]]);
		}
	}

	//check cell2Cell, to see some rules:
	//rule1 one cell's neighbor cells should not contain itself.
	//rule2 0<= one cell's neighbor cells < nTotalCell
	//rule3 one cell's neighbor cells should be different.
	//It should be noted that for O type grid, break 3 may be broken.
	void CheckCell2Cell(const int nTotalCell, const int * cell_number_of_each_cell, const int * offsetCell2Cell, const int * cell2cell){
		printf("Check cell2cell...\n");
		for (int cellID = 0; cellID < nTotalCell; cellID++){
			int offset = offsetCell2Cell[cellID];
			int numCellCells = cell_number_of_each_cell[cellID];
			for (int iCell = 0; iCell < numCellCells; iCell++){
				int cellCellID = cell2cell[offset + iCell];
				//check rule 1
				if (cellID == cellCellID) printf("rule 1 breaks in cellID %d in %d\n", cellID, iCell);
				//check rule 2
				if ((cellCellID < 0)||(cellCellID >= nTotalCell)) printf("rule 2 breaks in cellID %d in %d for %d\n", cellID, iCell, cellCellID);
			}

			//check rule 3
			int numComp = numCellCells;
			//At least two cells in cell2cell	
			while (numComp > 1){
				for (int jCell = 0; jCell < numComp-1; jCell++){
					int cellIDLoop = cell2cell[offset+jCell];
					int cellIDComp = cell2cell[offset+numComp-1];
					if (cellIDLoop == cellIDComp) printf("rule 3 breaks in cellID %d in %d and %d are both %d\n", cellID, jCell, numComp-1, cellIDLoop);
				}
				numComp--;
			}
		}

	}

	void ReOrderLevelCur(const int numSortCur, const int * cellIDSortCur, int & numVertex, int * NewCell){
		for (int iSort = 0; iSort < numSortCur; iSort++){
			int cellID = cellIDSortCur[iSort];
			//Avoid the same cell index in cellIDSortCur on the condition of O type grid
			//if (ReOrderCell[cellID] == -1){
			if (NewCell[cellID] == -1){
				numVertex++;
				NewCell[cellID] = numVertex;
			}
		}

	}

	void GetNextLevelForReOrder(const int & numSortCur, const int * cellIDSortCur, const int * cell_number_of_each_cell, const int * cell2cell, const int * offsetCell2Cell, int * ReOrderCell, int & numSortNext, int * cellIDSortNext){
		numSortNext = 0;
		for (int iCell = 0; iCell < numSortCur; iCell++){
			int cellCellID = cellIDSortCur[iCell];
			int numCellCells = cell_number_of_each_cell[cellCellID];
			int offset = offsetCell2Cell[cellCellID];
			for (int iCellCell = 0; iCellCell < numCellCells; iCellCell++){
				int cellID = cell2cell[offset + iCellCell];
				if (ReOrderCell[cellID] == -1){
					cellIDSortNext[numSortNext] = cellID;
					numSortNext++;
				}
			}
		}
	}
	//check rules for cell index
	//rule1 0 <= NewCell < nTotalCell
	//rule2 cell index should be different
	void CheckNewCell(const int nTotalCell, const int * NewCell){
		printf("Check NewCell ...\n");
		int * numCells = new int[nTotalCell];
		for (int iCell = 0; iCell < nTotalCell; iCell++){
			int cellID = NewCell[iCell];
			if ((cellID < 0)||(cellID > nTotalCell)) printf("CheckNewCell: cellID %d breaks rule1 on %d of NewCell\n", cellID, iCell);
			int remain = cellID % nTotalCell;
			numCells[remain]++;
		}
	
		for (int iCell = 0; iCell < nTotalCell; iCell++) {
			if (numCells[iCell] != 1) printf("CheckNewCell, cellID %d in NewCell breaks rule 2, has %d cells\n", iCell, numCells[iCell]);
		}	

		delete []numCells;
	}
	//compute face_number_of_each_cell and nTotalCellFaces
	//It should be noted that boundary faces are not considered here.
	void SetFaceNumberOfEachCell(const int nTotalFace, const int nTotalCell, const int * left_cell_of_face, const int * right_cell_of_face, int * face_number_of_each_cell, int & nTotalCellFaces){
		for (int iface = 0; iface < nTotalFace; iface++){
			//for left_cell_of_face
			int cellLabelLeft = left_cell_of_face[iface];
			if ((cellLabelLeft < 0)||(cellLabelLeft >= nTotalCell)) {
				printf("Error: cellLabelLeft = %d\n", cellLabelLeft);
				exit(1);
			}
			int cellLabelRight = right_cell_of_face[iface];
			//Boundary faces are not considered. In fact, on boundary face, cellLabelLeft is cell itself. cellLabelRight is -1, which is not necessary stored in cell2cell.
			if (cellLabelRight == -1) continue; //boundy face's right neighbor

			if ((cellLabelRight < 0)||(cellLabelRight >= nTotalCell)) {
				printf("Error: cellLabelRight = %d\n", cellLabelRight);
				exit(1);
			}
			face_number_of_each_cell[cellLabelLeft]++;
			nTotalCellFaces++;
			face_number_of_each_cell[cellLabelRight]++;
			nTotalCellFaces++;
		}		
	}
	//set offsetCell2Face by face_number_of_each_cell
	void SetOffSetCell2Face(const int nTotalCell, const int * face_number_of_each_cell, int * offsetCell2Face){
		offsetCell2Face[0] = 0;
		for (int icell = 1; icell < nTotalCell; icell++){
			offsetCell2Face[icell] = offsetCell2Face[icell-1] + face_number_of_each_cell[icell-1];
		}
	}
	//set cell2face by left_cell_of_face, right_cell_of_face, offsetCell2Face
	//It should be noted that boundary faces are not considered here.
	void SetCell2Face(const int nTotalFace, const int nTotalCell, const int * left_cell_of_face, const int * right_cell_of_face, int * offsetCell2Face, int * cell2face){
		for (int iface = 0; iface < nTotalFace; iface++){
			//for left_cell_of_face
			int cellLabelLeft = left_cell_of_face[iface];
			if ((cellLabelLeft < 0)||(cellLabelLeft >= nTotalCell)) {
				printf("Error: cellLabelLeft = %d\n", cellLabelLeft);
				exit(1);
			}

			int cellLabelRight = right_cell_of_face[iface];
			//boundary face is not considered.
			if (cellLabelRight == -1) continue; //boundy face's right neighbor

			if ((cellLabelRight < 0)||(cellLabelRight >= nTotalCell)) {
				printf("Error: cellLabelRight = %d\n", cellLabelRight);
				exit(1);
			}
			int cell2faceOffset = offsetCell2Face[cellLabelLeft];
			cell2face[cell2faceOffset] = iface;
			//offsetCell2Face update for the next set of face
			offsetCell2Face[cellLabelLeft]++;
			cell2faceOffset = offsetCell2Face[cellLabelRight];
			cell2face[cell2faceOffset] = iface;
			offsetCell2Face[cellLabelRight]++;
		}		
	}

	//set cell2Cell
	void SetCell2Cell(const int nTotalCell, const int * face_number_of_each_cell, const int * offsetCell2Face, const int * cell2face, const int * left_cell_of_face, const int * right_cell_of_face, int * cell2cell){
		for (int icell = 0; icell < nTotalCell; icell++){
			int sizeCellFace = face_number_of_each_cell[icell];
			int offsetCellFace = offsetCell2Face[icell];
			for (int iface = 0; iface < sizeCellFace; iface++){
				int faceID = cell2face[offsetCellFace + iface];
				int cellLabel = left_cell_of_face[faceID];
				if ((cellLabel < 0)||(cellLabel >= nTotalCell)) {
					printf("Error: cellLabel = %d\n", cellLabel);
					exit(1);
				}

				if (cellLabel != icell) cell2cell[offsetCellFace + iface] = cellLabel;

				cellLabel = right_cell_of_face[faceID];

				if (cellLabel == -1) continue; //boundy face's right neighbor
				if ((cellLabel < 0)||(cellLabel >= nTotalCell)) {
					printf("Error: cellLabel = %d\n", cellLabel);
					exit(1);
				}

				if (cellLabel != icell) cell2cell[offsetCellFace + iface] = cellLabel;
				
			}
		}
	}

	//compute bandwidth of each cell
	void computeBandwidth(const int nTotalCell, const int * face_number_of_each_cell, const int * offsetCell2Face, const int * cell2cell){
		int * bandwidth_of_each_cell = new int[nTotalCell];		
		for (int icell = 0; icell < nTotalCell; icell++){
			int sizeCellFace = face_number_of_each_cell[icell];
			int offsetCellFace = offsetCell2Face[icell];

			//The max band between icell and neighbor cells
			int maxLocalBand = 0;
			for (int iface = 0; iface < sizeCellFace; iface++){
				int cellID = cell2cell[offsetCellFace+iface];
				int bandDiff = abs(icell - cellID);
				if (bandDiff > maxLocalBand) maxLocalBand = bandDiff;
				//test begin
					//if (icell == 0) printf("offsetCellFace = %d, iface = %d, cellID = %d, bandDiff = %d\n", offsetCellFace, iface, cellID, bandDiff);
				//test end
			}
			bandwidth_of_each_cell[icell] = maxLocalBand;
		}
		//statistics of bandwidth_of_each_cell
		int maxGlobalBand = 0; //maximum bandwidth
		int maxCellID = 0; //cellID of maximum bandwidth
		/*
		int bandDist[33]; //distribution of bandwidth from 1 to 32 and more
		for (int i = 0; i < 33; i++){
			bandDist[i] = 0;
		}

		for (int icell = 0; icell < nTotalCell; icell++){
			int bandCell = bandwidth_of_each_cell[icell];
			if (bandCell <= 0) {
				printf("Error: bandwidth<=0 on cell %d, bandwidth=%d\n", icell, bandCell);
				exit(1);
			}
			if (bandCell > maxGlobalBand) {
				maxGlobalBand = bandCell;
				maxCellID = icell;
			}
			if (bandCell > 32) bandDist[32]++;
			else bandDist[bandCell-1]++;
		}
		
		printf("maxGlobalBand = %d on cellID %d\n", maxGlobalBand, maxCellID);
		int nTotalDist = 0;
		for (int i = 0; i < 33; i++){
			nTotalDist += bandDist[i];
			printf("bandDist %d: %d in nTotalCells\n", i+1, bandDist[i]);
		}
		*/
		int Level = 20;
		//int gapLevel = nTotalCell/Level;
		int gapLevel = 6000/Level;
		printf("gapLevel = %d\n", gapLevel);
		int bandDist[20]; //distribution of bandwidth from 1 to 32 and more
		for (int i = 0; i < Level; i++){
			bandDist[i] = 0;
		}

		for (int icell = 0; icell < nTotalCell; icell++){
			int bandCell = bandwidth_of_each_cell[icell];
			if (bandCell <= 0) {
				printf("Error: bandwidth<=0 on cell %d, bandwidth=%d\n", icell, bandCell);
				exit(1);
			}
			if (bandCell > maxGlobalBand) {
				maxGlobalBand = bandCell;
				maxCellID = icell;
			}
			int ratio = floor(bandCell/gapLevel);
			if (ratio < 0) {
				printf("Error: ratio = %d, icell = %d, bandCell = %d, gapLevel=%d\n", ratio, icell, bandCell, gapLevel);
				exit(0);
			}
			if (ratio > (Level-1)) bandDist[Level-1]++;
			else bandDist[ratio]++;
		}
		
		printf("maxGlobalBand = %d on cellID %d\n", maxGlobalBand, maxCellID);
		int nTotalDist = 0;
		for (int i = 0; i < Level; i++){
			nTotalDist += bandDist[i];
			printf("bandDist %d: %d in nTotalCells\n", i, bandDist[i]);
		}
		printf("nTotalDist = %d, nTotalCell = %d\n", nTotalDist, nTotalCell);
		
		delete []bandwidth_of_each_cell;		
	}

	//Update face_number_of_each_face by swap method
	//Intput: nTotalCell, face_number_of_each_face, reflectCell
	//Output:SwapFace_number_of_each_cell
	//The function is used to provide a swap method to update a cell related varible
	//after cell reorder. In the future, cell2node can be update by the method
	void SwapFaceNumberOfEachFace(const int nTotalCell, const int * face_number_of_each_cell, const int * reflectCell, const int * NewCell, const int * NewFace_number_of_each_cell){
		//try swap method
		int * SwapFace_number_of_each_cell = new int[nTotalCell];
		//It should be notes that OrgCell and SwapCell are updated during swap
		//Finally, OrgCell should be equal to reflectCell and SwapCell should be equal to NewCell. 
		int * OrgCell = new int[nTotalCell];//Contianer index and old cellID reflection
		int * SwapCell = new int[nTotalCell]; //Old cellID and container index reflection
		for (int cellID = 0; cellID < nTotalCell; cellID++){
			SwapFace_number_of_each_cell[cellID] = face_number_of_each_cell[cellID];
			//Before swap, OrgCell is equal to SwapCell
			OrgCell[cellID] = cellID; 
			SwapCell[cellID] = cellID;
		}
/*
		//set SwapFace_number_of_each_cell by swap method with NewCell
		//The method fails due to inifinte loop by while 
		for (int cellID = 0; cellID < nTotalCell; cellID++){
			//just for test
				printf("For loop, cellID = %d\n", cellID);
				int level = 0;
			//test end
			//if (OrgCell[cellID] == NewCell[cellID]) continue;
			if (OrgCell[cellID] == reflectCell[cellID]) continue;
			//It should be noted that on the condition of OrgCell[cellID] != NewCell[cellID], it must be obtained that OrgCell[cellID] = cellID
			//int newCellID = NewCell[cellID];
			int newCellID = reflectCell[cellID];
			int currentID = OrgCell[cellID];
			//just for test
				printf("level %d, newCellID = %d\n", level, newCellID);
			//test end
			while (currentID != newCellID){
				//swap face number  of currentID and newCellID in SwapFace_number_of_each_cell
				int temp = SwapFace_number_of_each_cell[currentID];
				SwapFace_number_of_each_cell[currentID] = SwapFace_number_of_each_cell[newCellID];
				SwapFace_number_of_each_cell[newCellID] = temp;
				//swap OrgCell
				OrgCell[currentID] = newCellID;
				OrgCell[newCellID] = currentID;
				
				//reset currentID by newCellID
				//currentID = NewCell[currentID];
				currentID = newCellID;
				//reset newCellID
				//newCellID = NewCell[currentID];
				newCellID = reflectCell[currentID];
				//just for test
					level++;
					printf("cellID = %d, level %d, currentID = %d, newCellID = %d\n", cellID, level, currentID, newCellID);
				//test end
				
			}
		}
*/
		//After swap, OrgCell should be equal to reflectCell
		//After swap, SwapCell should be equal to NewCell
		for (int cellID = 0; cellID < nTotalCell; cellID++){
			//In one container, the corresponding old cellID is the same
			//in OrgCell and reflectCell. It means that no swap is required.
			if (OrgCell[cellID] == reflectCell[cellID]) continue;
			//targetID is the target of OrgCell[cellID]
			//SwapID is the container index corresponds to targetID
			//currentID is the real value of OrgCell[cellID]
			//Obviously, OrgCell[cellID] and OrgCell[SwapID] should be swap.
			//After swap, SwapCell should be updated. SwapCell[targetID] 
			//should be cellID. and SwapCell[currentID] should be SwapID.
			//In fact, Swap[targerID] and Swap[currentID] is also swap.
			int targetID = reflectCell[cellID];
			int SwapID = SwapCell[targetID];
			int currentID = OrgCell[cellID];
			if (SwapCell[currentID] != cellID) {
				//It means that something wrong in  swap of OrgCell and SwapCell
				printf("Error: Old cellID %d is not in the containerID %d\n", currentID, cellID);
				printf("Old cellID %d is in the containerID %d\n", currentID, SwapCell[currentID]);
				printf("Container %d owns old cellID %d\n", cellID, OrgCell[cellID]);
				exit(1);
			}
			//swap face_number_of_each_cell of cellID and SwapID
			int temp = SwapFace_number_of_each_cell[SwapID];
			SwapFace_number_of_each_cell[SwapID] = SwapFace_number_of_each_cell[cellID];
			SwapFace_number_of_each_cell[cellID] = temp;
			//swap OrgCell of cellID and SwapID
			OrgCell[cellID] = targetID;
			OrgCell[SwapID] = currentID;
			//swap SwapID of targerID and currentID
			SwapCell[targetID] = cellID;
			SwapCell[currentID] = SwapID;
		}
		//check results
		CheckSwapFaceNumber(nTotalCell, OrgCell, SwapCell, reflectCell, NewCell, SwapFace_number_of_each_cell, NewFace_number_of_each_cell);

		
	}

	//check results
	void CheckSwapFaceNumber(const int nTotalCell, const int * OrgCell, const int * SwapCell, const int * reflectCell, const int * NewCell, const int * SwapFace_number_of_each_cell, const int * NewFace_number_of_each_cell){
		printf("Checking swap method for update face_number_of_each_cell just by itself...\n");
		for (int cellID = 0; cellID < nTotalCell; cellID++){
			//if (ReOrderCell[cellID] == -1) {
			//if (OrgCell[cellID] != NewCell[cellID]) {
			if (OrgCell[cellID] != reflectCell[cellID]) {
				printf("Error: cellID %d is not set yet\n", cellID);
				exit(1);
			}
			if (SwapFace_number_of_each_cell[cellID]!= NewFace_number_of_each_cell[cellID]) {
				printf("Error: on cellID %d, swap faces %d is not equal  to new faces %d\n", cellID, SwapFace_number_of_each_cell[cellID], NewFace_number_of_each_cell[cellID]);
				exit(1);
			}
		}
	}
	//Reorder cell label from nBoundFace to nTotalFace just with the order of faces in cell2face
	void cellReorderByCell2FaceOrder(const int nTotalCell, const int nTotalFace, const int * offsetCell2Face, const int * face_number_of_each_cell, const int * cell2face, const int nBoundFace, int labelFace, int * NewFaces){
		for (int cellID = 0; cellID < nTotalCell; cellID++){
			int offset = offsetCell2Face[cellID];
			int numFaces = face_number_of_each_cell[cellID];
			for (int iFace = 0; iFace < numFaces; iFace++){
				int faceID = cell2face[offset+iFace];
				//do not consider boundary faces
				if (faceID < nBoundFace) {
					printf("Error: faceID %d on cellID %d, iFace %d is boundary face\n", faceID, cellID, iFace);
					exit(1);
				}
				if (faceID >= nBoundFace){
					//faceID does not have label yet.
					if (NewFaces[faceID] == -1) {
						labelFace++;
						NewFaces[faceID] = labelFace;
					}
				}
			}
		}

		if (labelFace != nTotalFace-1) {
			printf("Error: labelFace %d is not equal to nTotalFace-1 %d\n", labelFace, nTotalFace-1);
			exit(1);
		}
	}
	//metric for computing left_cell_of_face and right_cell_of_face for GPU 
	void MetricLeftRightCellOfFace(const int nTotalFace, const int nBoundFace, const int * left_cell_of_face, const int * right_cell_of_face){
		int * diffCount = new int[129]; //0~127, 128 for larger
		//left part
		//int segAddress = 0;
		//start from nBoundFace, boundary faces are not considered
		int segAddress = nBoundFace;
		for (int i = 0; i < 129; i++){
			diffCount[i] = 0;
		}

		while (segAddress < nTotalFace){
				int localMax = left_cell_of_face[segAddress];
				int localMin = left_cell_of_face[segAddress];
				for (int offset = 0; offset < 32; offset++){
					int le = left_cell_of_face[segAddress];
					//printf("segAddress = %d, le = %d\n", segAddress, le);
					if (le < localMin) localMin = le;
					if (le > localMax) localMax = le;
					segAddress++;
				}
				int diff = localMax - localMin;
				if (diff > 128) diffCount[128]++;
				else diffCount[diff]++;
				//printf("localMax = %d, localMin = %d, diff = %d\n", localMax, localMin, diff);

			}
			//compute metric elements
			int totalNum = 0;
			for (int i = 0; i < 129; i++){
				totalNum += diffCount[i];
				//printf(" diffCount[%d] = %d\n", i, diffCount[i]);
			}
			for (int i = 0; i < 129; i++){
				double percent = double(diffCount[i])/totalNum;
				printf("i = %d, percent = %f\n", i, percent);
			}
			//printf("totalNum = %d, nTotalFace = %d\n", totalNum, nTotalFace);

			//right part
			//segAddress = 0;
			segAddress = nBoundFace;
			for (int i = 0; i < 129; i++){
				diffCount[i] = 0;
			}
			while (segAddress < nTotalFace){
				int localMax = right_cell_of_face[segAddress];
				int localMin = right_cell_of_face[segAddress];
				for (int offset = 0; offset < 32; offset++){
					int re = right_cell_of_face[segAddress];
					//printf("segAddress = %d, le = %d\n", segAddress, le);
					if (re < localMin) localMin = re;
					if (re > localMax) localMax = re;
					segAddress++;
				}
				int diff = localMax - localMin;
				if (diff > 128) diffCount[128]++;
				else diffCount[diff]++;
				//printf("localMax = %d, localMin = %d, diff = %d\n", localMax, localMin, diff);

			}
			//test
			totalNum = 0;
			for (int i = 0; i < 129; i++){
				totalNum += diffCount[i];
				//printf(" diffCount[%d] = %d\n", i, diffCount[i]);
			}
			for (int i = 0; i < 129; i++){
				double percent = double(diffCount[i])/totalNum;
				printf("i = %d, percent = %f\n", i, percent);
			}
			delete []diffCount;
		
	}
	//a caller for MetricLeftRightCellOfFace
	void CallMetricLeftRightCellOfFace(Base_Grid_Conn * gconn){
		printf("Program is running in CallMetricLeftRightCellOfFace\n");
		cgsize_t nTotalFace = gconn->GetNTotalFace();
		cgsize_t nBoundFace = gconn->GetNBoundFace();
		vector<cgsize_t> &left_cell_of_face  = *gconn->GetLeftCellOfFace();
    		vector<cgsize_t> &right_cell_of_face = *gconn->GetRightCellOfFace();
		int * Oldleft_cell_of_face = new int[nTotalFace];
		int * Oldright_cell_of_face = new int[nTotalFace];
		for (int iFace = 0; iFace < nTotalFace; iFace++){
			Oldleft_cell_of_face[iFace] = left_cell_of_face[iFace];
			Oldright_cell_of_face[iFace] = right_cell_of_face[iFace];
		}

		//call MetricLeftRightCellOfFace
		MetricLeftRightCellOfFace(nTotalFace, nBoundFace, Oldleft_cell_of_face, Oldright_cell_of_face);

	}

}
