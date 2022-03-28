// MemPool.h

#ifndef MEMPOOL_H
#define MEMPOOL_H
#include <memory>
#include <iostream>

template<typename C1, const int maxSize = 64, const int numPerNode = 1>
class MemHeap
{
protected:
	int freeIndex;
	C1 *heap;
	int *flag;
public:
  MemHeap() {
		heap = new C1[maxSize * numPerNode];
		flag = new int[maxSize];
		for (int i = 0; i < maxSize; i++) {
			flag[i] = i;
		}
	}
  MemHeap(const MemHeap &cpy) : freeIndex(cpy.freeIndex) {
		heap = new C1[maxSize * numPerNode];
		flag = new int[maxSize];
		for (int i = 0; i < maxSize; i++) {
			flag[i] = cpy.flag[i];
			for (int j = 0; j < numPerNode; j++) {
				heap[j + numPerNode * i] = cpy.heap[j + numPerNode * i];
			}
		}
	}
	~MemHeap() { delete[] heap; delete[] flag; }
public:
	C1 *get(int ptrId) {
		if(ptrId >= maxSize) return nullptr;
		return heap+numPerNode*ptrId; 
	}
	C1 *GetNew(int &ptrId) { 
		if(freeIndex == maxSize) {
			ptrId = maxSize;
			return nullptr;
		}
		else {
			ptrId = flag[freeIndex];
			freeIndex++;
			return heap+numPerNode * ptrId;
		}
	}
	bool del(const int &ptrId) { 
		for(int i = 0; i < freeIndex; i++) {
			if(flag[i] == ptrId) {
				freeIndex--;
				flag[i] = flag[freeIndex];
				flag[freeIndex] = ptrId;
				return true;
			}
		}
		return false;
	}
	const int &GetSize() const { return freeIndex; }
};

template<typename C1, const int maxSize, const int numPerNode = 1>
class MemPool
{
private:
	std::vector<std::shared_ptr<MemHeap<C1, maxSize, numPerNode>>> pool;
	int heapId;
public:
  MemPool() : heapId(0) {
		pool.push_back(std::shared_ptr<MemHeap<C1, maxSize, numPerNode>>(new MemHeap<C1, maxSize, numPerNode>));
	}
  MemPool(const MemPool &cpy) : heapId(cpy.heapId) {
    for(const auto &elem : cpy.pool) {
		  pool.push_back(std::shared_ptr<MemHeap<C1, maxSize, numPerNode>>(new MemHeap<C1, maxSize, numPerNode>(*elem)));
    }
	}
	C1* get(const int &dataId) const {
		int i1 = dataId/maxSize;
    if(i1 >= (int)pool.size() || pool[heapId] == nullptr) return nullptr;
		auto heap = pool[i1];
		return heap->get(dataId - i1 * maxSize);
	}
	C1 *GetNew(int &ptrId) {
		C1 *ptr = nullptr;
		for(; heapId < (int)pool.size() && pool[heapId] != nullptr; heapId++) {
			ptr = pool[heapId]->GetNew(ptrId);
			if(ptr != nullptr) {
				ptrId += maxSize * heapId;
				return ptr;
			}
		}
		pool.push_back(std::shared_ptr<MemHeap<C1, maxSize, numPerNode>>(new MemHeap<C1, maxSize, numPerNode>));
		ptr = pool[heapId]->GetNew(ptrId);
		ptrId += maxSize * heapId;
		return ptr;
	}
	bool del(int dataId) {
		int tempId = dataId / maxSize;
		if(tempId >= int(pool.size())) return false;
		auto heap = pool[tempId];
		if(heapId > tempId) heapId = tempId;
		return heap->del(dataId - tempId * maxSize);
	}
  void clear() { 
    pool.resize(1);
    heapId = 0;
    pool[heapId] = std::shared_ptr<MemHeap<C1, maxSize, numPerNode>>(new MemHeap<C1, maxSize,numPerNode>);
  }
};

template<typename C1, const int buffSize>
class MemGen
{
private:
  MemGen(const MemGen&);
private:
  const size_t numPerNode;
  int num;
  C1 *ptrBuffer[buffSize];
public:
  MemGen(const size_t &npn) : numPerNode(npn), num(0) {}
  ~MemGen() {
    //cout << "~MemGenTT() is running." << endl;
    for(int i = 0; i < num; i++) delete[] ptrBuffer[i];
  }
  void recycle(C1 *ptr) {
    if(ptr == nullptr) return;
    if(num == buffSize) {
      delete[] ptr;
      return;
    }
    ptrBuffer[num] = ptr;
    num++;
  }
  C1 *GetNew() { 
    if(num == 0) return new C1[numPerNode];
    num--;
    return ptrBuffer[num];
  }
  const size_t GetNumPerNode() const { return numPerNode; }
};

#endif