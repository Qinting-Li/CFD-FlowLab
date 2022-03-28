/////////////////////////////////////////////////
// author: Chen, Songze.  email: jacksongze@hust.edu.cn
/////////////////////////////////////////////////
// collaborator 1 :
/////////////////////////////////////////////////
// test 1 :
///////////////////////////////////////////////

// Range.h
#ifndef RANGE_H
#define RANGE_H

#include <set>
#include <map>
#include <vector>
#include <list>
#include <fstream>
#include <cstring>
#include "MemPool.h"
#include "BinBase.h"

#ifdef _DEBUG
#define CAL_ID(range, ...) CalIDSaveMacro(range, ##__VA_ARGS__) 
#else 
#define CAL_ID(range, ...) CalIDMacro(range, ##__VA_ARGS__)
#endif

inline size_t CalIDMacro(const int *range, const int *loc, const int &dim) {
  int ans = 0;
  int factor = 1;
  for(int i = dim - 1; i >= 0; i--) {
    factor = range[i + dim] - range[i];
    ans = loc[i] - range[i] + factor*ans;
  }
  return ans;
}

namespace nnwlbm{

	template<const size_t dim>
	class Range : public nnwlbm::DataClass
	{
	public:
		virtual std::string name() const { return "Range" + std::to_string(dim); } 
	protected:
		int range[2 * dim];
	public:
		inline const int &left(const size_t &i) const { return range[i]; }
		inline const int &right(const size_t &i) const { return range[i + dim]; }
		inline int &left(const size_t &i) { return range[i]; }
		inline int &right(const size_t &i) { return range[i + dim]; }
		inline int length(const size_t &i) const { return right(i) - left(i); }
		static Range CalIntersection(const Range &r1, const Range &r2) {
			Range ans;
			for (int i = 0; i < dim; i++) {
				ans.left(i) = (r1.left(i) > r2.left(i)) ? r1.left(i) : r2.left(i);
				ans.right(i) = (r1.right(i) > r2.right(i)) ? r2.right(i) : r1.right(i);
			}
			return std::move(ans);
		}
		static Range CalExtent(const Range &r1, const Range &r2) {
			Range ans;
			for (int i = 0; i < dim; i++) {
				ans.left(i) = (r1.left(i) > r2.left(i)) ? r2.left(i) : r1.left(i);
				ans.right(i) = (r1.right(i) > r2.right(i)) ? r1.right(i) : r2.right(i);
			}
			return ans;
		}
		static Range merge(const Range &r1, const Range &r2) {
			int counter = 0;
			int loc = 0;
			for (int i = 0; i < dim; i++) {
				if (r1.left(i) != r2.left(i) || r1.right(i) != r2.right(i)) {
					counter++;
					loc = i;
				}
			}
			if (counter > 1) return Range();
			if ((r1.left(loc) - r2.right(loc)) * (r1.right(loc) - r2.left(loc)) <= 0) return CalExtent(r1, r2);
			else return Range();
		}
		static Range enlarge(Range r1, const int &bufferSize) {
			for (int i = 0; i < dim; i++) {
				r1.left(i) -= bufferSize;
				r1.right(i) += bufferSize;
			}
			return r1;
		}
		// remove the overlapping part from r1.
		// return the remaining as a range list.
		static std::list<Range> substract(const Range &r1, const Range &r2) {
			std::list<Range> ans;
			Range ruler = CalIntersection(r1, r2);
			if (ruler.check()) {
				int length[dim], loc[dim];
				for (int i = 0; i < dim; i++) {
					loc[i] = i;
					length[i] = ruler.right(i) - ruler.left(i);
				}
				for (int i = 0; i < dim; i++) {
					// find min direction
					for (int j = i + 1; j < dim; j++) {
						if (length[j] < length[i]) {
							std::swap(loc[j], loc[i]);
							std::swap(length[j], length[i]);
						}
					}
					Range tempRangeL = ruler;
					Range tempRangeR = ruler;
					tempRangeL.left(loc[i]) = r1.left(loc[i]);
					tempRangeL.right(loc[i]) = ruler.left(loc[i]);
					if (tempRangeL.check()) ans.push_front(tempRangeL);
					tempRangeR.left(loc[i]) = ruler.right(loc[i]);
					tempRangeR.right(loc[i]) = r1.right(loc[i]);
					if (tempRangeR.check()) ans.push_front(tempRangeR);
					ruler.left(loc[i]) = r1.left(loc[i]);
					ruler.right(loc[i]) = r1.right(loc[i]);
				}
			}
			else {
				ans.push_front(r1);
			}
			return ans;
		}
		// mode == 0, split the larger range among r1 and r2
		// mode == 1, split range r1
		// mode == 2, split range r2
		// return the summation as a range list.
		static std::list<Range> sum(const Range &r1, const Range &r2, size_t mode = 0) {
			std::list<Range> ans;
			Range ruler = CalIntersection(r1, r2);
			if (!ruler.check()) {
				if (r1.check()) ans.push_front(r1);
				if (r2.check()) ans.push_front(r2);
			}
			else if (r1 == ruler) {
				ans.push_front(r2);
			}
			else if (r2 == ruler) {
				ans.push_front(r1);
			}
			else {
				int size1 = r1.GetSize();
				int size2 = r2.GetSize();
				Range splitRange;
				switch (mode) {
				case 1:
					ans.push_front(r2);
					splitRange = r1;
				case 2:
					ans.push_front(r1);
					splitRange = r2;
				default:
					if (size1 > size2) {
						ans.push_front(r2);
						splitRange = r1;
					}
					else {
						ans.push_front(r1);
						splitRange = r2;
					}
				};
				int length[dim], loc[dim];
				for (int i = 0; i < dim; i++) {
					loc[i] = i;
					length[i] = ruler.right(i) - ruler.left(i);
				}
				for (int i = 0; i < dim; i++) {
					// find min direction
					for (int j = i + 1; j < dim; j++) {
						if (length[j] < length[i]) {
							std::swap(loc[j], loc[i]);
							std::swap(length[j], length[i]);
						}
					}
					Range tempRangeL = ruler;
					Range tempRangeR = ruler;
					tempRangeL.left(loc[i]) = splitRange.left(loc[i]);
					tempRangeL.right(loc[i]) = ruler.left(loc[i]);
					if (tempRangeL.check()) ans.push_front(tempRangeL);
					tempRangeR.left(loc[i]) = ruler.right(loc[i]);
					tempRangeR.right(loc[i]) = splitRange.right(loc[i]);
					if (tempRangeR.check()) ans.push_front(tempRangeR);
					ruler.left(loc[i]) = splitRange.left(loc[i]);
					ruler.right(loc[i]) = splitRange.right(loc[i]);
				}
			}
			return ans;
		}
		static std::list<Range> gen_range_buffer_zoom(const Range &r1, const int &bufferSize) {
			if (bufferSize > 0) return substract(enlarge(r1, bufferSize), r1);
			else return std::list<Range>();
		}
	public:
		Range() { for (auto &elem : range) elem = 0; }
		Range(const int *data) {
			memcpy(range, data, 2 * dim * sizeof(int));
		}
		Range(int *data) {
			memcpy(range, data, 2 * dim * sizeof(int));
		}
		Range(std::vector<int> &data) {
			if (data.size() < 2 * dim) {
				for (int i = 0; i < data.size(); i++) range[i] = data[i];
				for (int i = data.size(); i < 2 * dim; i++) range[i] = 0;
			}
			else memcpy(range, data.data(), 2 * dim * sizeof(int));
		}
		Range(const std::vector<int> &data) {
			if (data.size() < 2 * dim) {
				for (int i = 0; i < data.size(); i++) range[i] = data[i];
				for (int i = data.size(); i < 2 * dim; i++) range[i] = 0;
			}
			else memcpy(range, data.data(), 2 * dim * sizeof(int));
		}
		Range(const Range &cpy) {
			memcpy(range, cpy.range, 2 * dim * sizeof(int));
		}
		template<typename ... T>
		Range(const T & ... arg) {
			int test[] = { arg... };
			int bound = 2 * dim < sizeof(test) ? 2 * dim : sizeof(test);
			for (int i = 0; i < bound; i++) { range[i] = test[i]; }
		}
		~Range() {}

		bool check() const {
			for (int i = 0; i < dim; i++) {
				if (left(i) >= right(i)) return false;
			}
			return true;
		}
		const int GetSize() const {
			int ans = 1;
			for (int i = 0; i < dim; i++) {
				ans *= right(i) - left(i);
			}
			return ans;
		}
		constexpr static inline size_t GetDim() { return dim; }

		int CalID(const int *loc) const {
			//return CalIDMacro(range, loc, dim);
			return CAL_ID((*this), loc);
		}
		bool CalIJK(int *loc, size_t id) const {
			for (size_t i = 0; i < GetDim(); i++) {
				const int length = (right(i) - left(i));
				loc[i] = id % length + left(i);
				id /= length;
			}
			const size_t lastI = GetDim() - 1;
			if (loc[lastI] < left(lastI) || loc[lastI] >= right(lastI)) return false;
			else return true;
		}
		int IniCarry(int *loc) const {
			if (!check()) return 1;
			for (int i = 0; i < GetDim(); i++) {
				loc[i] = left(i);
			}
			return 0;
		}
		int AddOne(int *loc) const {
			int i = 0;
			for (; i < dim; i++) {
				loc[i] += 1;
				if (loc[i] == right(i)) {
					loc[i] = left(i);
				}
				else break;
			}
			return (i < dim) ? 0 : 1;
		}
		int AddOne(int *loc, const int iDim) const {
			int i = 0;
			for (; i < dim; i++) {
				int ii = (i + iDim) % dim;
				loc[ii] += 1;
				if (loc[ii] == right(ii)) {
					loc[ii] = left(ii);
				}
				else break;
			}
			return (i < dim) ? 0 : 1;
		}
		//inline const int* data() const { return range; }
		//inline const int& operator[](const int& i) const { return range[i]; } 
		//inline int& operator[](const int& i) { return range[i]; }
		Range &operator = (const Range &right) {
			memcpy(this->range, right.range, 2 * dim * sizeof(int));
			return *this;
		}
		bool operator==(const Range &right) const {
			for (int i = 0; i < 2 * dim; i++) {
				if (range[i] != right.range[i]) return false;
			}
			return true;
		}
		bool operator!=(const Range &right) const { return !(*this == right); }
		bool operator<(const Range &right) const {
			for (int i = 0; i < 2 * dim; i++) {
				if (range[i] < right.range[i]) return true;
				else if (range[i] == right.range[i]) continue;
				else return false;
			}
			return false;
		}
		bool operator>(const Range &right) const {
			return !(*this == right || *this < right);
		}
		bool IsNeighbor(Range r2, const int &bufferSize = 1) const {
			Range intersection1 = CalIntersection(*this, r2);
			r2.AddBuffer(bufferSize);
			Range intersection2 = CalIntersection(*this, r2);
			return (!(intersection1.check())) && (intersection2.check());
		}
		bool IsContaining(const Range &r2) const {
			Range intersection = CalIntersection(*this, r2);
			return r2 == intersection;
		}
		bool IsContaining(const int *loc) const {
		  for (size_t i = 0; i < dim; i++) {
			  if (loc[i] < left(i) || loc[i] >= right(i)) return false;
		  }
		  return true;;
		}
		void AddBuffer(int bufferSize) {
			for (int i = 0; i < dim; i++) {
				left(i) -= bufferSize;
				right(i) += bufferSize;
			}
		}
		// set range
		template<class ... T>
		void SetLeft(const T& ... arg) {
			int test[] = { arg... };
			int bound = dim < sizeof(test) ? dim : sizeof(test);
			for (int i = 0; i < bound; i++) { left(i) = test[i]; }
		}
		template<class ... T>
		void SetRight(const T& ... arg) {
			int test[] = { arg... };
			int bound = dim < sizeof(test) ? dim : sizeof(test);
			for (int i = 0; i < bound; i++) { right(i) = test[i]; }
		}
		inline void GetLeft(int *s) const {
			for (size_t i = 0; i < GetDim(); i++) {
				s[i] = left(i);
			}
		}
		inline void GetRight(int *e) const {
			for (size_t i = 0; i < GetDim(); i++) {
				e[i] = right(i);
			}
		}
		void SetStart(const int* s) {
			for (size_t i = 0; i < dim; i++) { left(i) = s[i]; }
		}
		void SetEnd(const int* e) {
			for (size_t i = 0; i < dim; i++) { right(i) = e[i] + 1; }
		}
	private:
		virtual bool OutputBinPrivate(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const // 必须
		{
			const size_t temp = GetDim();
			ostr.write((const char*)(range), 2 * temp * sizeof(int));
			return true;
		}

    friend std::ostream &operator << (std::ostream &out, const Range &rg) {
      out << "[(";
      out << rg.left(0);
      for(size_t i = 1; i < dim; i++) {
        out << ' ' << rg.left(i);
      }
      out << "), (";
      out << rg.right(0);
      for(size_t i = 1; i < dim; i++) {
        out << ' ' << rg.right(i);
      }
      out << ")]";
      return out;
    }
	};

  //std::ostream& operator << (std::ostream& out, const Range<2>& rg) {
  //  out << "[( ";
  //  for(size_t i = 0; i < 2; i++) {
  //    out << rg.left(i) << " ";
  //  }
  //  out << "), (";
  //  for(size_t i = 0; i < 2; i++) {
  //    out << rg.right(i) << " ";
  //  }
  //  out << ")]";
  //  return out;
  //}

	template<int dim, int h = 1>
	Range<dim - h> DegenerateTail(const Range<dim>&rr, int loc = 0) {
		Range<dim - h> ans;
		for (int i = 0; i < dim - h; i++) {
			ans.left(i) = rr.left(i);
			ans.right(i) = rr.right(i);
		}
	}
	template<int dim, int h = 1>
	Range<dim - h> DegenerateHead(const Range<dim>&rr) {
		Range<dim - h> ans;
		for (int i = 0; i < dim - h; i++) {
			ans.left(i) = rr.left(i + h);
			ans.right(i) = rr.right(i + h);
		}
	}

	template<const size_t dim>
	size_t CalIDSaveMacro(const Range<dim>& rr, const int* loc) {
		int ans = 0;
		for (size_t i = dim - 1; i >= 0; i--) {
			if (loc[i] < rr.left(i) || loc[i] >= rr.right(i)) {
        std::cout << "CalIDSaveMacro() : out of range" << std::endl;
#ifdef _WIN32
        system("pause");
#endif
        exit(1);
      }
			ans = loc[i] - rr.left(i) + rr.length(i)*ans;
      if(i == 0) break;
		}
		return ans;
	}

    //typename std::enable_if<sizeof...(arg) != dim, void>::type
	template<size_t dim, class ... T>
  size_t  CalIDSaveMacro(
    const Range<dim> &rr, 
    const T& ... arg
  )
  {
		int loc[] = { arg... };
    if(sizeof...(arg) != dim) {
      std::cout << "CalIDSaveMacro() : sizeof...(arg) != dim" << std::endl;
#ifdef _WIN32
      system("pause");
#endif
      exit(1);
    }
		return CalIDSaveMacro(rr, loc);
	}

	inline size_t CalIDMacro(const Range<1> &rr, const int *loc) {
		return loc[0] - rr.left(0);
	}
	inline size_t CalIDMacro(const Range<2> &rr, const int *loc) {
		return loc[0] - rr.left(0) + rr.length(0)*(loc[1] - rr.left(1));
	}
	inline size_t CalIDMacro(const Range<3> &rr, const int *loc) {
		return loc[0] - rr.left(0) + (rr.right(0) - rr.left(0))*
			(loc[1] - rr.left(1) + (rr.right(1) - rr.left(1))*
			(loc[2] - rr.left(2)));
	}
	inline size_t CalIDMacro(const Range<4> &rr, const int *loc) {
		return loc[0] - rr.left(0) + (rr.right(0) - rr.left(0))*
			(loc[1] - rr.left(1) + (rr.right(1) - rr.left(1))*
			(loc[2] - rr.left(2) + (rr.right(2) - rr.left(2))*
				(loc[3] - rr.left(3))));
	}

	template<int dim, class ... T>
	inline size_t CalIDMacro(const Range<dim> &rr, const T& ... arg) {
		int loc[] = { arg... };
		return CalIDMacro(rr, loc);
	}

	// Constructor: the input range is nodeRange.
	template<const size_t dim>
	class RangeIndex
	{
	private:
		Range<dim> cellRange;
		Range<dim> nodeRange;
	public:
		enum { inputCellRange, inputNodeRange };
		RangeIndex(const Range<dim> &r, const int &mode /*= inputCellRange*/)
			: cellRange(r), nodeRange(r)
		{
			for (size_t i = 0; i < GetDim(); i++) {
				nodeRange.right(i) += 1;
			}
			switch (mode) {
			case inputCellRange:
				break;
			case inputNodeRange:
			  for (size_t i = 0; i < GetDim(); i++) {
					cellRange.right(i) -= 1;
					nodeRange.right(i) -= 1;
				}
				break;
			default:
				std::cout << "RangeIndex : wrong mode." << std::endl;
			}
		}
		inline size_t GetDim() const { return cellRange.GetDim(); }
		inline const Range<dim> &GetCellRange() const {
			return cellRange;
		}
		inline const Range<dim> &GetNodeRange() const {
			return nodeRange;
		}
		inline std::vector<Range<dim>> GetFaceRange() const {
			std::vector<Range<dim>> ans;
			for (int i = 0; i < GetDim(); i++) {
				Range<dim> fr(cellRange);
				fr.right(i) += 1;//在相应维的方向上+1
				ans.push_back(fr);
			}
			return ans;
		}
		inline size_t CalFaceID(const int *loc, const size_t &direction) const {
			return direction + (int)GetDim()*CAL_ID(nodeRange, loc);
		}
		inline size_t CalNodeID(const int *loc) const {
			return CAL_ID(nodeRange, loc);
		}
		inline size_t CalCellID(const int *loc) const {
			return CAL_ID(cellRange, loc);
		}
		inline bool CalFaceIJK(int *loc, size_t &direction, const size_t &id) const {
			direction = id % GetDim();
			return nodeRange.CalIJK(loc, id / GetDim());
		}
		inline bool CalNodeIJK(int *loc, const size_t &id) const {
			return nodeRange.CalIJK(loc, id);
		}
		bool CalCellIJK(int *loc, const size_t &id) const {
			return cellRange.CalIJK(loc, id);
		}
	};

	// non-overlapping range set
	template<const size_t dim>
	class RangeSet {
	protected:
		Range<dim> extent;
		std::set<Range<dim>> rangeSet;
		void UpdateExtent() {
			if (rangeSet.empty()) {
				extent = Range<dim>();
				return;
			}
			extent = *(rangeSet.begin());
			for (const auto &elem : rangeSet) {
				for (int i = 0; i < dim; i++) {
					extent.left(i) = (extent.left(i) > elem.left(i)) ? elem.left(i) : extent.left(i);
					extent.right(i) = (extent.right(i) > elem.right(i)) ? extent.right(i) : elem.right(i);
				}
			}
		}
	public:
		RangeSet() {}
		bool check() const {
			for (const auto &elem1 : rangeSet) {
				for (const auto &elem2 : rangeSet) {
					if (elem1 == elem2) continue;
					if (CalIntersection(elem1, elem2).check()) return false;
				}
			}
			return true;
		}
		void add(const Range<dim> &range) {
			std::list<Range<dim>> rangeList;
			rangeList.push_back(range);
			if (!extent.check()) extent = range;
			else extent = Range<dim>::CalExtent(extent, range);
			for (const Range<dim> &elem : rangeSet) {
				for (auto iter = rangeList.begin(); iter != rangeList.end(); ) {
					bool isSplitted = Range<dim>::CalIntersection(*iter, elem).check();
					if (isSplitted) {
						auto subRangeList = Range<dim>::substract(*iter, elem);
						rangeList.splice(rangeList.end(), subRangeList);
						rangeList.erase(iter);
						iter = rangeList.begin();
					}
					else iter++;
				}
			}
			for (const auto &tempR : rangeList) {
				rangeSet.insert(tempR);
			}
		}
		void substract(const Range<dim> &range) {
			std::list<Range<dim>> delList;
			std::list<Range<dim>> addList;
			for (const Range<dim> &elem : rangeSet) {
				bool isSplitted = Range<dim>::CalIntersection(range, elem).check();
				if (isSplitted) {
					auto subRangeList = Range<dim>::substract(elem, range);
					addList.splice(addList.end(), subRangeList);
					delList.push_back(elem);
				}
			}
			for (const Range<dim> &elem : delList) {
				rangeSet.erase(elem);
			}
			auto iter = rangeSet.begin();
			for (const Range<dim> &elem : addList) {
				iter = rangeSet.insert(iter, elem);
			}
			UpdateExtent();
		}
		const std::set<int, Range<dim>> &GetRangeSet() const { return rangeSet; }
	};

	template<typename C, const size_t dim>
	class MaskRange : public nnwlbm::DataClass{
	protected:
		Range<dim> range;
		std::vector<C> msk;
	public:
		template<class ... T>
		MaskRange(const T& ... arg) : range(arg...), msk(range.GetSize(), 0) {  }
		explicit MaskRange(const std::vector<int> &data, const int &m = 0) : range(data), msk(range.GetSize(), m) {  }
		explicit MaskRange(const Range<dim> &rg, const int &m = 0) : range(rg), msk(rg.GetSize(), m) {  }
		explicit MaskRange(const MaskRange &cpy) : range(cpy.range), msk(cpy.msk) {  }
		MaskRange(const Range<dim> &rg, const std::vector<C> &mskIn) : range(rg), msk(mskIn) {  }

		const std::vector<C> &mask() const { return msk; }

		virtual void SetMask(const size_t &id, const C &val) { msk[id] = val; }
		C GetMask(const size_t &id) const { return msk[id]; }
		template <class Compare>
		size_t count_num(Compare comp) const {
			size_t i = 0;
			for (const C &elem : msk) if (comp(elem)) i++;
			return i;
		}
		Range<dim> GetRange() const { return range; }
		MaskRange<int, dim> GetSubMaskRange(const Range<dim> &rg) const
		{
			MaskRange<int, dim> ans(rg);
			int loc[dim];
			for (int carry = rg.IniCarry(loc); carry == 0; carry = rg.AddOne(loc)) {
				int id = CAL_ID(range, loc);
				int nId = CAL_ID(rg, loc);
				ans.SetMask(nId, this->GetMask(id));
			}
			return ans;
		}
	public:
		virtual std::string name() const { return "MaskRange" + std::to_string(sizeof(C)) +std::to_string(dim); } // 必须
	private:
		virtual bool OutputBinPrivate(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const // 必须
		{
			int temprange[2 * dim];
			for (int i = 0; i < dim; i++) {
				temprange[i] = range.left(i);
				temprange[i + dim] = range.right(i);
			}
			ostr.write((const char*)(temprange), 2 * dim * sizeof(int));
			size_t temp = msk.size(); 
			ostr.write((const char*)(&temp), 1 * sizeof(size_t)); 
			ostr.write((const char*)(msk.data()), msk.size() * sizeof(C)); 
			return true;
		}
	};

	// MaskData : the index only represents the range. The variable is returned as a pointer.
	struct Color {
		bool operator() (const int &i) const { return i == 1 || i == 0; }
	};

  template<typename T, typename BinaryComparator>
  class Comparator {
  public:
    Comparator(T secd) : second(secd), binaryComparator(BinaryComparator()) {}
    bool operator() (const T& first) const { return binaryComparator(first, second); }
  private:
    T second;
    BinaryComparator binaryComparator;
  }; 

	struct ColorInterface {
		virtual int GetColor(const size_t &id) const = 0;
		virtual void SetColor(const size_t &id, const int &val) = 0;
	};

	template<typename C, const size_t dim, typename Compare = Color, const int ndef = -99999>
	class MaskData : public MaskRange<int, dim>, public ColorInterface {
	private:
		MemGen<C, 64> *memGen;//连续开64个内存，放置指针
	protected:
		std::vector<C*> var;
		Compare compare;
	public:
		using MaskRange<int, dim>::msk;
		using MaskRange<int, dim>::range;
		using MaskRange<int, dim>::GetRange;
		MaskData(const size_t &vn, const MaskRange<int, dim> &mask, Compare comp = Compare())
			: MaskRange<int, dim>(mask), memGen(new MemGen<C, 64>(vn)), compare(comp) {//vn表示每一个节点上对应的数据个数
			var.resize(range.GetSize(), nullptr);
			for (int i = range.GetSize() - 1; i >= 0; i--) {
				if (comp(msk[i])) var[i] = memGen->GetNew();
			}
		}
		MaskData(const size_t &vn, const Range<dim> &rg, Compare comp)
			: MaskRange<int, dim>(rg), memGen(new MemGen<C, 64>(vn)), compare(comp) {
			var.resize(range.GetSize(), nullptr);
			for (int i = range.GetSize() - 1; i >= 0; i--) {
				SetMask(i, ndef);
			}
		}
		MaskData(const MaskData<C, dim, Compare> &md)
			: MaskRange<int, dim>(md.GetRange()), compare(md.compare),
			memGen(new MemGen<C, 64>(md.GetVarNum()))
		{
			var.resize(range.GetSize(), nullptr);
			BeFilledWith(md);
		}
		MaskData(MaskData<C, dim, Compare>&& md)
			: MaskRange<int, dim>(md.GetRange()), compare(md.compare),
			memGen(md.memGen)
		{
			md.memGen = nullptr;
			var = std::move(md.var);
		}

		~MaskData() {
			for (C *elem : var) memGen->recycle(elem);
			delete memGen;
		}
		operator const MaskRange<int, dim>&() const { return *this; }

		const Compare& get_comp() const { return compare; }
		const size_t GetVarNum() const { return memGen->GetNumPerNode(); }
		bool IsValid(const size_t &id) const { return compare(GetColor(id)); }//判断id是否是需要的
    RangeIndex<dim> GetNodeIndex() const { return RangeIndex<dim>(GetRange(), RangeIndex<dim>::inputNodeRange); }
		virtual void SetColor(const size_t &id, const int &val) {
			SetMask(id, val);
		}
		virtual int GetColor(const size_t &id) const {
			if (id < 0) return ndef;
			return MaskRange<int, dim>::GetMask(id);
		}
		virtual void SetMask(const size_t& id, const int &val) {
			msk[id] = val;
			if (compare(msk[id])) {
				if (var[id] == nullptr) var[id] = memGen->GetNew();
			}
			else {
				if (var[id] != nullptr) {
					memGen->recycle(var[id]);
					var[id] = nullptr;
				}
			}
		}
		C *GetVar(const size_t &id) { return var[id]; }
		const C *GetVar(const size_t &id) const { return var[id]; }
		// mode == 0, transfer the ownership
		void transfer_data_from(MaskData<C, dim, Compare> &md) {
			Range<dim> intersection = Range<dim>::CalIntersection(this->GetRange(), md.GetRange());
			int loc[dim];
			for (int carry = intersection.IniCarry(loc); carry == 0; carry = intersection.AddOne(loc))
			{
				int thisId = CAL_ID(GetRange(), loc);
				int mdId = CAL_ID(md.GetRange(), loc);
				msk[thisId] = md.msk[mdId];
				memGen->recycle(var[thisId]);
				var[thisId] = md.var[mdId];
				md.msk[mdId] = ndef;
				md.var[mdId] = nullptr;
			}
		}
		// mode == 1, cpy the data，var of md dont change
		void BeFilledWith(const MaskData<C, dim, Compare> &md) {
			Range<dim> intersection = Range<dim>::CalIntersection(GetRange(), md.GetRange());
			int loc[dim];
			for (int carry = intersection.IniCarry(loc); carry == 0; carry = intersection.AddOne(loc))
			{
				size_t thisId = CAL_ID(GetRange(), loc);
				size_t mdId = CAL_ID(md.GetRange(), loc);
				SetMask(thisId, md.GetMask(mdId));
				if (var[thisId] != nullptr && md.var[mdId] != nullptr) {
					for (int i = 0; i < GetVarNum(); i++) {
						var[thisId][i] = md.var[mdId][i];
					}
				}
			}
		}
	};

	template<typename C, const size_t dim, class Compare = Color>
	class MaskDataMap : public RangeSet<dim> {
	private:
		std::map<Range<dim>, std::shared_ptr<MaskData<C, dim, Compare>>> mskDataMap;
		const int varNum;
	public:
		using RangeSet<dim>::extent;
		using RangeSet<dim>::rangeSet;
		using RangeSet<dim>::UpdateExtent;
		MaskDataMap(const int& vn) : varNum(vn) {}
		void add(const Range<dim> &range) {
			std::list<Range<dim>> rangeList;
			rangeList.push_back(range);
			if (!extent.check()) extent = range;
			else extent = Range<dim>::CalExtent(extent, range);

			for (const Range<dim> &elem : rangeSet) {
				for (auto iter = rangeList.begin(); iter != rangeList.end();) {
					bool isSplitted = Range<dim>::CalIntersection(*iter, elem).check();
					if (isSplitted) {
						auto subRangeList = Range<dim>::substract(*iter, elem);
						rangeList.splice(rangeList.end(), subRangeList);
						rangeList.erase(iter);
						iter = rangeList.begin();
					}
					else iter++;
				}
			}
			for (const auto &tempR : rangeList) {
				rangeSet.insert(tempR);
				mskDataMap.emplace(tempR, std::make_shared<MaskData<C, dim, Compare>>(varNum, tempR, Compare()));
			}
		}
		void add(MaskData<C, dim, Compare> &md) {
			Range<dim> range = md.GetRange();
			std::list<Range<dim>> rangeList;
			rangeList.push_back(range);
			if (!extent.check()) extent = range;
			else extent = Range<dim>::CalExtent(extent, range);

			for (const Range<dim> &elem : rangeSet) {
				for (auto iter = rangeList.begin(); iter != rangeList.end(); iter++) {
					bool isSplitted = Range<dim>::CalIntersection(*iter, elem).check();
					if (isSplitted) {
						auto subRangeList = Range<dim>::substract(*iter, elem);
						rangeList.splice(rangeList.end(), subRangeList);
						rangeList.erase(iter);
						iter = rangeList.begin();
					}
					else iter++;
				}
			}
			for (const auto &tempR : rangeList) {
				rangeSet.insert(tempR);
				auto piece = std::make_shared<MaskData<C, dim, Compare>>(varNum, tempR, Compare());
				piece->transfer_data_from(md);
				mskDataMap.emplace(tempR, piece);
			}
		}

		void substract(const Range<dim> &range) {
			std::list<Range<dim>> delList;
			std::list<std::shared_ptr<MaskData<C, dim, Compare>>> addList;
			for (const Range<dim> &elem : rangeSet) {
				bool isSplitted = Range<dim>::CalIntersection(range, elem).check();
				if (isSplitted) {
					auto subRangeList = Range<dim>::substract(elem, range);
					auto originalPiece = mskDataMap.at(elem);
					for (const auto& tempR : subRangeList) {
						auto piece = std::make_shared<MaskData<C, dim, Compare>>(varNum, tempR, Compare());
						piece->transfer_data_from(*originalPiece);
						addList.push_back(piece);
					}
					delList.push_back(elem);
				}
			}
			for (const Range<dim> &elem : delList) {
				rangeSet.erase(elem);
				mskDataMap.erase(elem);
			}
			for (auto& piece : addList) {
				rangeSet.insert(piece->GetRange());
				mskDataMap.emplace(piece->GetRange(), piece);
			}
			UpdateExtent();
		}

		const MaskData<C, dim, Compare> &GetMaskData(const Range<dim> &rg) const { return *(mskDataMap.at(rg)); }
		MaskData<C, dim, Compare> &GetMaskData(const Range<dim> &rg) { return *(mskDataMap.at(rg)); }
	};

	template<const size_t dim, class Compare = Color>
	std::shared_ptr<MaskRange<int, dim>> BuildDiscreteLevelSet(const MaskRange<int, dim> &mr, Compare comp) {
		Range<dim> rg = mr.GetRange();
		int loc[dim];
		std::vector<int> ls(rg.GetSize(), -1);
		std::set<size_t> pending;
		std::set<size_t> completed;
		for (int carry = rg.IniCarry(loc); carry == 0; carry = rg.AddOne(loc))
		{
			size_t id = CAL_ID(rg, loc);
			int msk = mr.GetMask(id);
			auto spread = [&]() {
				const size_t locId = *(pending.begin());
				pending.erase(pending.begin());
				int loc1[dim];
				rg.CalIJK(loc1, locId);
				Range<dim> neighborRange;
				neighborRange.SetStart(loc1);
				neighborRange.SetEnd(loc1);
				neighborRange.AddBuffer(1);
				neighborRange = Range<dim>::CalIntersection(neighborRange, rg);
				for (int c1 = neighborRange.IniCarry(loc1); c1 == 0; c1 = neighborRange.AddOne(loc1))
				{
					size_t neighborId = CAL_ID(rg, loc1);
					if (comp(mr.GetMask(neighborId))) {
						auto iter = completed.find(neighborId);
						if (iter == completed.end()) {
							ls[neighborId] = 0;
							pending.insert(neighborId);
							completed.insert(neighborId);
						}
						continue;
					}
					if (ls[neighborId] == -1 || ls[neighborId] > ls[locId] + 1) {
						ls[neighborId] = ls[locId] + 1;
						pending.insert(neighborId);
					}
				}
			};
			if (comp(msk)) {
				auto iter = completed.find(id);
				if (iter == completed.end()) {
					ls[id] = 0;
					pending.insert(id);
					completed.insert(id);
				}
				while (!pending.empty()) {
					spread();
				}
				break;
			}
		}
		return std::make_shared<MaskRange<int, dim>>(rg, ls);
	}

	template<const size_t dim, class Compare>
	Range<dim> get_subExtent(const MaskRange<int, dim> &mrg, Compare comp)
	{
		Range<dim> point;
		Range<dim> ans;
		int loc[dim];
		Range<dim> rg = mrg.GetRange();
		for (int carry = rg.IniCarry(loc); carry == 0; carry = rg.AddOne(loc)) {
			int id = CAL_ID(rg, loc);
			if (comp(mrg.GetMask(id))) {
				point.SetStart(loc); point.SetEnd(loc);
				ans = Range<dim>::CalExtent(point, ans);
			}
		}
		return ans;
	}

	template<const size_t dim>
	std::shared_ptr<MaskRange<int, dim>> GetSubMaskRange(const MaskRange<int, dim> &mr, const Range<dim> &rg)
	{
		std::shared_ptr<MaskRange<int, dim>> ans = std::make_shared<MaskRange<int, dim>>(rg);
		int loc[dim];
		auto mrRg = mr.GetRange();
		for (int carry = rg.IniCarry(loc); carry == 0; carry = rg.AddOne(loc)) {
			int id = CAL_ID(mrRg, loc);
			int nId = CAL_ID(rg, loc);
			ans->SetMask(nId, mr.GetMask(id));
		}
		return ans;
	}

	// class Comp1 is used to build level set
	template<const size_t dim, class Comp1>
	std::pair<std::list<std::shared_ptr<MaskRange<int, dim>>>, bool>
		SplitMaskRange(
			const MaskRange<int, dim> &mr, // the maskRange to be split
			Comp1 comp1, // to build level set
			int minSize // minimun ruler size
		) {
		auto ans = std::make_pair(std::list<std::shared_ptr<MaskRange<int, dim>>>(), true);
		auto levelSet = BuildDiscreteLevelSet(mr, comp1);
		MaskRange<int, dim> &tempMr = *levelSet;
		Range<dim> rg = tempMr.GetRange();
		int loc[dim], maxD, maxId;
		maxD = 0; maxId = 0;
		for (int carry = rg.IniCarry(loc); carry == 0; carry = rg.AddOne(loc))
		{
			int id = CAL_ID(rg, loc);
			if (maxD < tempMr.GetMask(id)) {
				maxD = tempMr.GetMask(id);
				maxId = id;
			}
		}
		if (minSize < maxD) {
			// if split
			std::set<int> completed;
			Range<dim> innerRg;
			// find an inner point whose neighbors possess most number of maxD 
			std::list<int> idList;
			idList.push_back(maxId);
			int preCount = -1, nId = maxId;
			int count = 0;
			while (preCount < count && idList.size() != 0) {
				preCount = count;
				count = 0;
				nId = idList.back();
				idList.pop_back();
				completed.insert(nId);
				rg.CalIJK(loc, nId);
				innerRg.SetStart(loc); innerRg.SetEnd(loc);
				innerRg.AddBuffer(1);
				innerRg = Range<dim>::CalIntersection(rg, innerRg);
				for (int c1 = innerRg.IniCarry(loc); c1 == 0; c1 = innerRg.AddOne(loc))
				{
					int id = CAL_ID(rg, loc);
					if (tempMr.GetMask(id) == maxD) {
						count++;
						if (completed.find(id) == completed.end()) {
							idList.push_front(id);
						}
					}
				}
			}
			// generate a rectangle area covered by points with maxD.
			rg.CalIJK(loc, nId);
			innerRg.SetStart(loc); innerRg.SetEnd(loc);
			for (int iDim = 0; iDim < dim; iDim++) {
				// extend to left
				bool isValid = true;
				while (isValid) {
					Range<dim> marginRg(innerRg);
					marginRg.left(iDim) = innerRg.left(iDim) - 1;
					marginRg.right(iDim) = marginRg.left(iDim) + 1;
					marginRg = Range<dim>::CalIntersection(rg, marginRg);
					if (!marginRg.check()) {
						isValid = false;
					}
					for (int c1 = marginRg.IniCarry(loc); c1 == 0; c1 = marginRg.AddOne(loc)) {
						int id = CAL_ID(rg, loc);
						if (tempMr.GetMask(id) != maxD) {
							isValid = false;
							break;
						}
					}
					if (isValid) {
						innerRg.left(iDim) -= 1;
					}
				}
				// extend to right
				isValid = true;
				while (isValid) {
					Range<dim> marginRg(innerRg);
					marginRg.left(iDim) = innerRg.right(iDim);
					marginRg.right(iDim) = marginRg.left(iDim) + 1;
					marginRg = Range<dim>::CalIntersection(rg, marginRg);
					if (!marginRg.check()) {
						isValid = false;
					}
					for (int c1 = marginRg.IniCarry(loc); c1 == 0; c1 = marginRg.AddOne(loc)) {
						int id = CAL_ID(rg, loc);
						if (tempMr.GetMask(id) != maxD) {
							isValid = false;
							break;
						}
					}
					if (isValid) {
						innerRg.right(iDim) += 1;
					}
				}
			}
			innerRg.AddBuffer(maxD - 1);
			auto rgList = Range<dim>::substract(rg, innerRg);
			// generate maskrange from rgList
			for (const auto& subRg : rgList) {
				ans.first.push_back(GetSubMaskRange(mr, subRg));
			}
			ans.second = true;
		}
		else {
			// if not split
			ans.first.push_back(GetSubMaskRange(mr, mr.GetRange()));
			ans.second = false;
		}
		return ans;
	}
	template<const size_t dim, class Comp1>
	std::list<std::shared_ptr<MaskRange<int, dim>>>
		partition_maskRange(
			const MaskRange<int, dim>& mr, // the maskRange to be split
			Comp1 comp1, // to build level set
			int minSize // minimun ruler size
		) {
		std::list<std::shared_ptr<MaskRange<int, dim>>> completed;
		std::list<std::shared_ptr<MaskRange<int, dim>>> pending;
		auto smr = SplitMaskRange(mr, comp1, minSize);
		if (!smr.second) {
			completed.splice(completed.begin(), smr.first);
		}
		else {
			pending.splice(pending.begin(), smr.first);
		}
		while (pending.size() != 0) {
			auto smr = SplitMaskRange(*(pending.back()), [](int v) { return v == 0; }, 0);
			pending.pop_back();
			if (!smr.second) {
				completed.splice(completed.begin(), smr.first);
			}
			else {
				pending.splice(pending.begin(), smr.first);
			}
			return completed;
		}
	}

	//Binary istream of range
	template<const size_t dim>
	class RangeRead : public nnwlbm::BinInput
	{
	public:
		virtual std::string name() const { return "Range" + std::to_string(dim); }
	private:
		virtual std::shared_ptr<nnwlbm::DataClass> InputBinPrivate(
			std::ifstream &istr, /**< [in,out]二进制文件输入流 */
			const std::string &filename, /**< [in]读入的文件名 */
			const std::string &path = ".\\" /**< [in]文件路径，默认路径为当前目录 */
		) const {
			int range[2 * dim];
			istr.read((char*)(range), 2 * dim * sizeof(int));
			//make_shared()函数可以接受最多10个参数，然后把它们传递给类型T的构造函数，创建一个shared_ptr<T>的对 象并返回。
			return std::make_shared<Range<dim>>(range);
		}
	};

	//Binary read of maskrange
	template<typename C, const size_t dim>
	class MaskRangeRead : public nnwlbm::BinInput
	{
	public:
		virtual std::string name() const { return "MaskRange" + std::to_string(sizeof(C)) + std::to_string(dim); }
	private:
		virtual std::shared_ptr<nnwlbm::DataClass> InputBinPrivate(
			std::ifstream &istr, /**< [in,out]二进制文件输入流 */
			const std::string &filename, /**< [in]读入的文件名 */
			const std::string &path = ".\\" /**< [in]文件路径，默认路径为当前目录 */
		) const {
			int range[2 * dim];
			istr.read((char*)(range), 2 * dim * sizeof(int));
			nnwlbm::Range<dim> temprange(range);
			// 临时内存空间
			size_t tempSize;
			// 读入int数据个数
			istr.read((char*)(&tempSize), sizeof(size_t));
			// 分配内存
			std::vector<C> CVec(tempSize);
			// 读入int数据
			istr.read((char*)(CVec.data()), tempSize * sizeof(C));
			// 相当于对MaskRange进行初始化
			return std::make_shared<MaskRange<C, dim>>(temprange, CVec);
		}
	};

  class LinkedData
  {
  public:
    virtual void *allocate(const size_t &id) = 0;
    friend class abd;
  };

	template<typename C, const size_t dim>
  class LinkedMaskRange
  {
  private:
    std::set<LinkedData*> dataSet;
		Range<dim> range;
		std::vector<C> msk;
  public:
		template<class ... T>
		LinkedMaskRange(const T& ... arg) 
      : range(arg...), msk(range.GetSize(), 0) {  }

		explicit LinkedMaskRange(
      const std::vector<int> &data, 
      const int &m = 0
    ) : range(data), msk(range.GetSize(), m) {  }

		explicit LinkedMaskRange(
      const Range<dim> &rg, 
      const int &m = 0
    ) : range(rg), msk(rg.GetSize(), m) {  }

		LinkedMaskRange(const LinkedMaskRange &cpy) 
      : range(cpy.range), msk(cpy.msk), dataSet(cpy.dataSet) {  }

		LinkedMaskRange(
      const Range<dim> &rg, 
      const std::vector<C> &mskIn
    ) : range(rg), msk(mskIn) {  }

		const std::vector<C> &mask() const { return msk; }

		inline void SetColor(const size_t &id, const C &val) { 
      msk[id] = val; 
      for(auto data : dataSet) {
        data->allocate(id);
      }
    }

		inline C GetColor(const size_t &id) const { return msk[id]; }

		template <class Compare>
		size_t count_num(Compare comp) const {
			size_t i = 0;
			for (const C &elem : msk) if (comp(elem)) i++;
			return i;
		}

		constexpr static const size_t GetDim() { return Range<dim>::GetDim(); }

		inline Range<dim> GetRange() const { return range; }

    inline int GetSize() const { return range.GetSize(); }

    inline void link(LinkedData *ptr) { 
      dataSet.emplace(ptr); 
    }
    inline void unlink(LinkedData *ptr) { 
      dataSet.erase(ptr); 
    }
  };

	template<typename C, const size_t dim, typename Compare = Color, const int ndef = -99999>
	class LinkedMaskData : public LinkedData {
	private:
		MemGen<C, 64> *memGen;//连续开64个内存，放置指针
    LinkedMaskRange<int, dim> &hostMask;
	protected:
		std::vector<C*> var;
		Compare compare;
	public:
		LinkedMaskData(const size_t &vn, LinkedMaskRange<int, dim> &mask, Compare comp = Compare())
			: hostMask(mask), memGen(new MemGen<C, 64>(vn)), compare(comp) {//vn表示每一个节点上对应的数据个数
      hostMask.link(this);
			var.resize(hostMask.GetSize(), nullptr);
			for (int i = (int)var.size() - 1; i >= 0; i--) {
				if (IsValid(i)) var[i] = memGen->GetNew();
			}
		}

		LinkedMaskData(const LinkedMaskData<C, dim, Compare> &md)
			: hostMask(md.hostMask), compare(md.compare), memGen(new MemGen<C, 64>(md.GetVarNum()))
		{
      hostMask.link(this);
			var.resize(hostMask.GetSize(), nullptr);
			BeFilledWith(md);
		}

		LinkedMaskData(LinkedMaskData<C, dim, Compare>&& md)
			: hostMask(md.hostMask), compare(md.compare), memGen(md.memGen)
		{
      hostMask.unline(&md);
      hostMask.link(this);
			md.memGen = nullptr;
			var = std::move(md.var);
		}

		~LinkedMaskData() {
			for (C *elem : var) memGen->recycle(elem);
			delete memGen;
		}

		const LinkedMaskRange<int, dim> &GetMask() const { return hostMask; }

		const Compare &get_comp() const { return compare; }

		const size_t GetVarNum() const { return memGen->GetNumPerNode(); }

		bool IsValid(const size_t &id) const { return compare(GetColor(id)); }//判断id是否是需要的

    Range<dim> GetRange() const { return hostMask.GetRange(); }
    
    RangeIndex<dim> GetNodeIndex() const { return RangeIndex<dim>(hostMask.GetRange(), RangeIndex<dim>::inputNodeRange); }

		void SetColor(const size_t &id, const int &val) {
			hostMask.SetColor(id, val);
		}

		int GetColor(const size_t &id) const {
			return hostMask.GetColor(id);
		}

    size_t count_num(Compare comp) const { return hostMask.count_num(comp); }

		inline C *GetVar(const size_t &id) { return var[id]; }
		inline const C *GetVar(const size_t &id) const { return var[id]; }

    void BeFilledWith(const LinkedMaskData<C, dim, Compare> &md) {
      Range<dim> intersection = Range<dim>::CalIntersection(GetRange(), md.GetRange());
      int loc[dim];
      for (int carry = intersection.IniCarry(loc); carry == 0; carry = intersection.AddOne(loc))
      {
        size_t thisId = CAL_ID(GetRange(), loc);
        size_t mdId = CAL_ID(md.GetRange(), loc);
        SetColor(thisId, md.GetColor(mdId));
        if (var[thisId] != nullptr && md.var[mdId] != nullptr) {
          for (int i = 0; i < GetVarNum(); i++) {
            var[thisId][i] = md.var[mdId][i];
          }
        }
      }
    }
  private:
		void* allocate(const size_t &id) {
			if (IsValid(id)) {
				if (var[id] == nullptr) var[id] = memGen->GetNew();
			}
			else {
				if (var[id] != nullptr) {
					memGen->recycle(var[id]);
					var[id] = nullptr;
				}
			}
      return (void*)(var[id]);
		}
	};

}
#endif

// end of IndexTool.h
