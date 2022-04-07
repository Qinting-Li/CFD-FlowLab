#include "D2Q9.h"
#include "D3Q15.h"
#include "TimeMarching.h"
#include "UniformGrid.h"
#include "LatticeMPI.h"
#include <iostream> 
#include <mpi.h> 

const uint32_t DirectionMask::directionInBinary[32] = { 
  0b1,                      // 0
  0b10,                     // 1
  0b100,                    // 2 
  0b1000,                   // 3 
  0b10000,                  // 4 
  0b100000,                 // 5 
  0b1000000,                // 6
  0b10000000,               // 7
  0b100000000,              // 8
  0b1000000000,             // 9
  0b10000000000,                      // 10
  0b100000000000,                     // 11
  0b1000000000000,                    // 12
  0b10000000000000,                   // 13
  0b100000000000000,                  // 14
  0b1000000000000000,                 // 15
  0b10000000000000000,                // 16
  0b100000000000000000,               // 17
  0b1000000000000000000,              // 18
  0b10000000000000000000,             // 19
  0b100000000000000000000,                  // 20
  0b1000000000000000000000,                  // 21
  0b10000000000000000000000,                  // 22
  0b100000000000000000000000,                  // 23
  0b1000000000000000000000000,                  // 24
  0b10000000000000000000000000,                  // 25
  0b100000000000000000000000000,                  // 26
  0b1000000000000000000000000000,                  // 27
  0b10000000000000000000000000000,                  // 28
  0b100000000000000000000000000000,                  // 29
  0b1000000000000000000000000000000,                  // 30
  0b10000000000000000000000000000000                  // 31
};

//class RequestHeap
RequestHeap::RequestHeap() {
  maxSize = 0;
  increment = 26;
  freeIndex = 0;
  completeIndex = 0;
  flag = NULL;
  status = NULL;
  heap = NULL;
}

void RequestHeap::reset() {
  maxSize = 0;
  increment = 26;
  freeIndex = 0;
  completeIndex = 0;
  if(flag != NULL) {
    delete[] flag;
    flag = NULL;
  }
  if(status != NULL) {
    delete[] status;
    status = NULL;
  }
  if(heap != NULL) {
    delete[] heap;
    heap = NULL;
  }
}

RequestHeap::~RequestHeap() {
  reset();
}

MPI_Request* RequestHeap::get_request(int i) {
  if(i < 0 || i >= freeIndex) return NULL;
  return heap+i; 
}

MPI_Request* RequestHeap::GetNew_request() {
  if(freeIndex + 2 >= maxSize) {
    //		if(!is_some_complete()) {
    if(!expand()) return NULL;
    //		}
  }
  freeIndex = freeIndex + 1;
  //	MPI_Request_free(heap+freeIndex-1);
  return heap+freeIndex - 1;
}

const MPI_Status* RequestHeap::get_complete_status() {
  if(completeIndex <= 0) return NULL;
  return status; 
}

int RequestHeap::is_some_complete() {
  if(freeIndex == 0) return completeIndex = 0;
  MPI_Testsome(freeIndex, heap, &completeIndex, flag, status);
  int i = 0;
  int j = 0;
  int counter = 0;

  for(j = 0; j < completeIndex; j++) {
    for(; i < flag[j]; i++) {
      heap[counter] = heap[i];
      counter++;
    }
    i++;
  }

  for(; i < freeIndex; i++) {
    heap[counter] = heap[i];
    counter++;
  }

  freeIndex = freeIndex - completeIndex;
  return completeIndex;
}

bool RequestHeap::test_all_complete() {
  //if(freeIndex == 0) return true;
  bool isAllCompleted;
  int boolFlag;
  MPI_Testall(freeIndex, heap, &boolFlag, status);
  isAllCompleted = (boolFlag != 0);
  if(isAllCompleted) {
    //for(int i = 0; i < freeIndex; i++) {
    //  if(*(heap+i) != MPI_REQUEST_NULL) {
    //    MPI_Request_free(heap+i);
    //  }
    //}
    completeIndex = freeIndex;
    freeIndex = 0;
  }
  return isAllCompleted;
}

void RequestHeap::WaitAllComplete(){
  using std::cout;
  using std::endl;
  if(freeIndex == 0) return;
  int flag = MPI_Waitall(freeIndex, heap, status);

  if(flag != MPI_SUCCESS) {
    cout << "MPI_Waitall returns flag " << flag << ". Look in status." << endl;
    for(int i = 0; i < freeIndex; i++) {
      cout << "i == " << i << "\n" << "heap[i] " << heap[i] << " status[i] " 
        << status[i].MPI_SOURCE <<' ' << status[i].MPI_TAG << ' ' << status[i].MPI_ERROR << endl;
    }
  }
  completeIndex = freeIndex;
  freeIndex = 0;
}

bool RequestHeap::expand() {
  MPI_Request* tempReq = new MPI_Request[maxSize + increment];
  bool* tempBool = new bool[maxSize + increment];
  int* tempInt = new int[maxSize + increment];
  MPI_Status* tempStatus = new MPI_Status[maxSize + increment];
  if((!tempReq) || (!tempBool) || (!tempInt) || (!tempStatus) ) return false;
  if(heap) {
    for(int i = 0; i < freeIndex; i++) {
      tempReq[i] = heap[i];
    }
    delete[] heap;
  }
  heap = tempReq;
  if(status) {
    for(int i = 0; i < completeIndex; i++) {
      tempStatus[i] = status[i];
    }
    delete[] status;
  }
  status = tempStatus;
  if(flag) delete[] flag;
  flag = tempInt;

  maxSize = maxSize + increment;
  return true;
}

void unpack_lattice_pdf(
  double* pdf,
  DirectionMask& dirMask,
  const void* buffer, 
  int bufferSize, 
  int* pos, 
  MPI_Comm comm,
  const int& maxVelNum = 32
) {
  MPI_Unpack(buffer, bufferSize, pos, dirMask.data(), 1, MPI_UINT32_T, comm); 
  for(size_t dPos = 0; dirMask.next_direction(dPos, maxVelNum); ) {
    MPI_Unpack(buffer, bufferSize, pos, pdf + dPos - 1, 1, MPI_DOUBLE, comm); 
  }
}

void pack_lattice_pdf(
  void* buffer, 
  int bufferSize, 
  int* pos, 
  const double* pdf,
  const DirectionMask& dirMask,
  MPI_Comm comm,
  const int& maxVelNum = 32
) {
  MPI_Pack(dirMask.data(), 1, MPI_UINT32_T, buffer, bufferSize, pos, comm); 
  for(size_t dPos = 0; dirMask.next_direction(dPos, maxVelNum); ) {
    MPI_Pack(pdf + dPos - 1, 1, MPI_DOUBLE, buffer, bufferSize, pos, comm); 
  }
}

