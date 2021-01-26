#ifdef USE_ARRAY_OPT
#include "Math_BasisFunction.h"
#include "DataStruct_Array.h"
#include "PHOpt.h"
#include <string.h>
#include "Constants.h"

#ifdef WIN32
//! Shield compiling warning C4068 (by Zhangjian).
//! This warning is because of some optimization pragma not fit for MSVC compiler.
#pragma warning(disable: 4068)
#endif
namespace PHSPACE {
//template void fmsub(PHArray<RFloat, 3> dst, const PHArray<RFloat, 3> &src1, const PHArray<RFloat, 3> &src2);

// static
//inline
void local_update_msubadd_xyz(int ist, int iend, RDouble*  base_dstx_1, RDouble *base_dstx1_1, RDouble * __restrict base_srcx1_1, RDouble * __restrict base_srcx2_1, RDouble*  base_dsty_1, RDouble *base_dsty1_1, RDouble * __restrict base_srcy1_1, RDouble*  base_dstz_1, RDouble *base_dstz1_1, RDouble * __restrict base_srcz1_1) {
	#pragma ivdep
	#pragma vector always	
#pragma unroll(8)
#pragma prefetch base_srcx1_1:1:256
#pragma prefetch base_srcy1_1:1:256
#pragma prefetch base_srcz1_1:1:256
#pragma prefetch base_dstx_1:1:256
#pragma prefetch base_dstx1_1:1:256
#pragma prefetch base_dsty_1:1:256
#pragma prefetch base_dsty1_1:1:256
#pragma prefetch base_dstz_1:1:256
#pragma prefetch base_dstz1_1:1:256
	for (int i = ist ; i < iend; i++) {
		RDouble base_src = base_srcx2_1[i];
		RDouble tmp_x = base_srcx1_1[i] * base_src;
		RDouble tmp_y = base_srcy1_1[i] * base_src;
		RDouble tmp_z = base_srcz1_1[i] * base_src;
		
		base_dstx_1[i] -= tmp_x;
		base_dstx1_1[i] += tmp_x;
		base_dsty_1[i] -= tmp_y;
		base_dsty1_1[i] += tmp_y;
		base_dstz_1[i] -= tmp_z;
		base_dstz1_1[i] += tmp_z;
	}
}

void local_XfnMultiAreaAddXfnMultiArea(int ist, int iend, RFloat * __restrict dstx, RFloat * __restrict dsty, RFloat * __restrict dstz, const RFloat *srcx1, const RFloat *srcy1, const RFloat *srcz1, const RFloat *share1, const RFloat *srcx2, const RFloat *srcy2, const RFloat *srcz2, const RFloat *share2) {
	
	int i;
	#pragma ivdep
	#pragma vector always	
#pragma unroll(8)
	for (i = ist; i < iend; i++) {
        RFloat share1_tmp = share1[i];
        RFloat share2_tmp = share2[i];
		dstx[i] = srcx1[i] * share1_tmp + srcx2[i] * share2_tmp;
		dsty[i] = srcy1[i] * share1_tmp + srcy2[i] * share2_tmp;
		dstz[i] = srcz1[i] * share1_tmp + srcz2[i] * share2_tmp;
	}
}
void local_XfnMultiAreaAddXfnMultiArea(int ist, int iend, RFloat * __restrict dstx, RFloat * __restrict dsty, RFloat * __restrict dstz, const RFloat *srcx1, const RFloat *srcy1, const RFloat *srcz1, const RFloat *srcx2, const RFloat *srcy2, const RFloat *srcz2) {
	
	int i;
	#pragma ivdep
	#pragma vector always	
	#pragma unroll(8)
	for (i = ist; i < iend; i++) {
		dstx[i] = srcx1[i] + srcx2[i];
		dsty[i] = srcy1[i] + srcy2[i];
		dstz[i] = srcz1[i] + srcz2[i];
	}
}
void XfnMultiAreaAddXfnMultiArea(PHArray<RDouble, 3> &dstx, PHArray<RDouble, 3> &dsty, PHArray<RDouble, 3> &dstz, const PHArray<RDouble, 3> &srcx1, const PHArray<RDouble, 3> &srcy1, const PHArray<RDouble, 3>& srcz1, const PHArray<RDouble, 3> &share1, const PHArray<RDouble, 3> &srcx2, const PHArray<RDouble, 3> &srcy2, const PHArray<RDouble, 3>& srcz2, const PHArray<RDouble, 3> &share2) {

	int j, k;
	//! Declare *dstx_base_data = dstx.data_.
    declare_data(dstx); declare_data(dsty); declare_data(dstz);
    declare_data(srcx1); declare_data(srcy1); declare_data(srcz1); declare_data(share1);
    declare_data(srcx2); declare_data(srcy2); declare_data(srcz2); declare_data(share2);

	int stride_2 = dstx.stride_[2];
	int stride_1 = dstx.stride_[1];

    k = dstx.storage_.base(2);
    int kend = k + dstx.length_[2];
    for (; k < kend; k++) {
		//! Solve the length of third dimention; declare *dstx_base_2 = dstx_base_data + k * stride_2.
		declare_basex(dstx, 2, data, k); declare_basex(dsty, 2, data, k); declare_basex(dstz, 2, data, k);
		declare_basex(srcx1, 2, data, k); declare_basex(srcy1, 2, data, k); declare_basex(srcz1, 2, data, k); declare_basex(share1, 2, data, k);
		declare_basex(srcx2, 2, data, k); declare_basex(srcy2, 2, data, k); declare_basex(srcz2, 2, data, k); declare_basex(share2, 2, data, k);

		j = dstx.storage_.base(1);
		int jend = j + dstx.length_[1];
		for (; j < jend; j++) {
			//! Solve the length of third dimention; declare *dstx_base_1 = dstx_base_2 + j * stride_1.
		    declare_basex(dstx, 1, 2, j); declare_basex(dsty, 1, 2, j); declare_basex(dstz, 1, 2, j);
		    declare_basex(srcx1, 1, 2, j); declare_basex(srcy1, 1, 2, j); declare_basex(srcz1, 1, 2, j); declare_basex(share1, 1, 2, j);
		    declare_basex(srcx2, 1, 2, j); declare_basex(srcy2, 1, 2, j); declare_basex(srcz2, 1, 2, j); declare_basex(share2, 1, 2, j);

			int ist = dstx.storage_.base(0);
			int iend = ist + dstx.length_[0];
			//! dstx_base_1.
			local_XfnMultiAreaAddXfnMultiArea(ist, iend,
		        declare_base_namex(dstx, 1), declare_base_namex(dsty, 1), declare_base_namex(dstz, 1),
		        declare_base_namex(srcx1, 1), declare_base_namex(srcy1, 1), declare_base_namex(srcz1, 1), declare_base_namex(share1, 1),
		        declare_base_namex(srcx2, 1), declare_base_namex(srcy2, 1), declare_base_namex(srcz2, 1), declare_base_namex(share2, 1));
		}
    }

}
void XfnMultiAreaAddXfnMultiArea(PHArray<RDouble, 3> &dstx, PHArray<RDouble, 3> &dsty, PHArray<RDouble, 3> &dstz, const PHArray<RDouble, 3> &srcx1, const PHArray<RDouble, 3> &srcy1, const PHArray<RDouble, 3>& srcz1, const PHArray<RDouble, 3> &srcx2, const PHArray<RDouble, 3> &srcy2, const PHArray<RDouble, 3>& srcz2) {

	int j, k;
	//! Declare *dstx_base_data = dstx.data_.
    declare_data(dstx); declare_data(dsty); declare_data(dstz);
    declare_data(srcx1); declare_data(srcy1); declare_data(srcz1);
    declare_data(srcx2); declare_data(srcy2); declare_data(srcz2); 

	int stride_2 = dstx.stride_[2];
	int stride_1 = dstx.stride_[1];

    k = dstx.storage_.base(2);
    int kend = k + dstx.length_[2];
    for (; k < kend; k++) {
		//! Solve the length of third dimention; declare *dstx_base_2 = dstx_base_data + k * stride_2.
		declare_basex(dstx, 2, data, k); declare_basex(dsty, 2, data, k); declare_basex(dstz, 2, data, k);
		declare_basex(srcx1, 2, data, k); declare_basex(srcy1, 2, data, k); declare_basex(srcz1, 2, data, k); 
		declare_basex(srcx2, 2, data, k); declare_basex(srcy2, 2, data, k); declare_basex(srcz2, 2, data, k); 

		j = dstx.storage_.base(1);
		int jend = j + dstx.length_[1];
		for (; j < jend; j++) {
			//! Solve the length of third dimention; declare *dstx_base_1 = dstx_base_2 + j * stride_1.
		    declare_basex(dstx, 1, 2, j); declare_basex(dsty, 1, 2, j); declare_basex(dstz, 1, 2, j);
		    declare_basex(srcx1, 1, 2, j); declare_basex(srcy1, 1, 2, j); declare_basex(srcz1, 1, 2, j); 
		    declare_basex(srcx2, 1, 2, j); declare_basex(srcy2, 1, 2, j); declare_basex(srcz2, 1, 2, j); 

			int ist = dstx.storage_.base(0);
			int iend = ist + dstx.length_[0];
			//! dstx_base_1.
			local_XfnMultiAreaAddXfnMultiArea(ist, iend,
		        declare_base_namex(dstx, 1), declare_base_namex(dsty, 1), declare_base_namex(dstz, 1),
		        declare_base_namex(srcx1, 1), declare_base_namex(srcy1, 1), declare_base_namex(srcz1, 1), 
		        declare_base_namex(srcx2, 1), declare_base_namex(srcy2, 1), declare_base_namex(srcz2, 1) );
		}
    }

}
void local_QaddDqladdDqr(int ist, int iend, RDouble c1, RDouble c2, RDouble * __restrict dst1, RFloat * __restrict dst2, const RFloat * src1, const RFloat * src2, const RFloat * src3) {
    int i;
	register RDouble rt1, rt2, rt3;
	#pragma ivdep
	#pragma vector always
	for (i = ist; i < iend; i++) {
		rt1 = src1[i];
		rt2 = src2[i];
		rt3 = src3[i];
		dst1[i] = rt1 + c1 * rt2 + c2 * rt3;
		dst2[i] = rt1 - c2 * rt2 - c1 * rt3;
	}
}

void QaddDqladdDqr(RDouble c1, RDouble c2, PHArray<RDouble, 4> dst1, PHArray<RDouble, 4> dst2, const PHArray<RDouble, 4> &src1, const PHArray<RDouble, 4> &src2, const PHArray<RDouble, 4> &src3) {
	int j, k, m;
    declare_data(dst1); declare_data(dst2); 
    declare_data(src1); declare_data(src2); declare_data(src3); 

	int stride_3 = dst1.stride_[3];
	int stride_2 = dst1.stride_[2];
	int stride_1 = dst1.stride_[1];

    m = dst1.storage_.base(3);
    int mend = m + dst1.length_[3];
	for (; m < mend; m++) {
		declare_basex(dst1, 3, data, m); declare_basex(dst2, 3, data, m);
		declare_basex(src1, 3, data, m); declare_basex(src2, 3, data, m); declare_basex(src3, 3, data, m);

    	k = dst1.storage_.base(2);
    	int kend = k + dst1.length_[2];
		for (; k < kend; k++) {
			declare_basex(dst1, 2, 3, k); declare_basex(dst2, 2, 3, k);
			declare_basex(src1, 2, 3, k); declare_basex(src2, 2, 3, k); declare_basex(src3, 2, 3, k);

			j = dst1.storage_.base(1);
			int jend = j + dst1.length_[1];
			for (; j < jend; j++) {
				declare_basex(dst1, 1, 2, j); declare_basex(dst2, 1, 2, j);
				declare_basex(src1, 1, 2, j); declare_basex(src2, 1, 2, j); declare_basex(src3, 1, 2, j);

				int ist = dst1.storage_.base(0);
				int iend = ist + dst1.length_[0];
				local_QaddDqladdDqr(ist, iend, c1, c2,
				    declare_base_namex(dst1, 1), declare_base_namex(dst2, 1),
				    declare_base_namex(src1, 1), declare_base_namex(src2, 1), declare_base_namex(src3, 1));
			}
		}
	}

}
void local_QlrEqualDq(int ist,int iend,RFloat * __restrict dst1, RFloat * __restrict dst2, const RFloat * src1, const RFloat * src2) {
    int i;
	#pragma ivdep
	#pragma vector always
	for (i = ist; i < iend; i++) {
		dst1[i] = src1[i];
		dst2[i] = src2[i];
	}
}
void QlrEqualDq(PHArray<RDouble, 4> dst1, PHArray<RDouble, 4> dst2, const PHArray<RDouble, 4> &src1, const PHArray<RDouble, 4> &src2)
{
	int j, k, m;
    declare_data(dst1); declare_data(dst2); 
    declare_data(src1); declare_data(src2); 

	int stride_3 = dst1.stride_[3];
	int stride_2 = dst1.stride_[2];
	int stride_1 = dst1.stride_[1];

    m = dst1.storage_.base(3);
    int mend = m + dst1.length_[3];
	for (; m < mend; m++) {
		declare_basex(dst1, 3, data, m); declare_basex(dst2, 3, data, m);
		declare_basex(src1, 3, data, m); declare_basex(src2, 3, data, m);

    	k = dst1.storage_.base(2);
    	int kend = k + dst1.length_[2];
		for (; k < kend; k++) {
			declare_basex(dst1, 2, 3, k); declare_basex(dst2, 2, 3, k);
			declare_basex(src1, 2, 3, k); declare_basex(src2, 2, 3, k);

			j = dst1.storage_.base(1);
			int jend = j + dst1.length_[1];
			for (; j < jend; j++) {
				declare_basex(dst1, 1, 2, j); declare_basex(dst2, 1, 2, j);
				declare_basex(src1, 1, 2, j); declare_basex(src2, 1, 2, j);

				int ist = dst1.storage_.base(0);
				int iend = ist + dst1.length_[0];
				local_QlrEqualDq(ist, iend,
				    declare_base_namex(dst1, 1), declare_base_namex(dst2, 1),
				    declare_base_namex(src1, 1), declare_base_namex(src2, 1));
			}
		}
	}

}

void local_FourthMultiQ(RDouble fourth,int ist,int iend,RFloat * __restrict dst1, RFloat * __restrict src1, RFloat *__restrict src2, RFloat * __restrict src3, RFloat *__restrict src4) {
    int i;
#pragma ivdep
#pragma vector always
#pragma unroll(4)
#pragma prefetch src1:1:256
#pragma prefetch src2:1:256
#pragma prefetch src3:1:256
#pragma prefetch src4:1:256
	for (i = ist; i < iend; i++) {
         dst1[i]= fourth * (src1[i] + src2[i] + src3[i] + src4[i]);
	}
}

void QsubQ1(PHArray<RDouble, 4> dst1, const PHArray<RDouble, 4> &src1, const PHArray<RDouble, 4> &src2)
{
	int j, k, m;
    declare_data(dst1); 
    declare_data(src1); declare_data(src2); 

	int stride_3 = dst1.stride_[3];
	int stride_2 = dst1.stride_[2];
	int stride_1 = dst1.stride_[1];

    m = dst1.storage_.base(3);
    int mend = m + dst1.length_[3];
	for (; m < mend; m++) {
		declare_basex(dst1, 3, data, m);
		declare_basex(src1, 3, data, m); declare_basex(src2, 3, data, m);

    	k = dst1.storage_.base(2);
    	int kend = k + dst1.length_[2];
		for (; k < kend; k++) {
			declare_basex(dst1, 2, 3, k); 
			declare_basex(src1, 2, 3, k); declare_basex(src2, 2, 3, k);

			j = dst1.storage_.base(1);
			int jend = j + dst1.length_[1];
			for (; j < jend; j++) {
				declare_basex(dst1, 1, 2, j); 
				declare_basex(src1, 1, 2, j); declare_basex(src2, 1, 2, j);

				int ist = dst1.storage_.base(0);
				int iend = ist + dst1.length_[0];
				local_QsubQ1(ist, iend,
				    declare_base_namex(dst1, 1),
				    declare_base_namex(src1, 1), declare_base_namex(src2, 1));
			}
		}
	}

}

void local_TheReciprocalOfVolAddVol(int ist,int iend,RFloat * __restrict dst1, const RFloat * src1,const RFloat * src2) {
    int i;
#pragma ivdep
#pragma vector always
	for (i = ist; i < iend; i++) {
		dst1[i] = 1.0 / (src1[i] + src2[i]);
	}
}

void local_DqdxMultiWorkqm(int ist,int iend,RFloat * __restrict dst1, RFloat * __restrict dst2,RFloat * __restrict dst3,const RFloat * src1) {
    int i;
	register RDouble rt;
#pragma ivdep
#pragma vector always
	for (i = ist; i < iend; i++) {
        rt = src1[i];
		dst1[i] *= rt;
		dst2[i] *= rt;
		dst3[i] *= rt;
	}
}

void WorksxyzMultiQ(PHArray<RFloat, 3> dstx, PHArray<RFloat, 3> dstx1, PHArray<RFloat, 3> dsty, PHArray<RFloat, 3> dsty1,PHArray<RFloat, 3> dstz, PHArray<RFloat, 3> dstz1,const PHArray<RFloat, 3> &srcx1,const PHArray<RFloat, 3> &srcy1,const PHArray<RFloat, 3> &srcz1, const PHArray<RFloat, 3> &src2) 
{
	RFloat * data_dstx, *data_dstx1, *data_srcx1, *data_src2;
	RFloat * data_dsty, *data_dsty1, *data_srcy1;
	RFloat * data_dstz, *data_dstz1, *data_srcz1;

	data_dstx = dstx.data_;
	data_dstx1 = dstx1.data_;
	data_dsty = dsty.data_;
	data_dsty1 = dsty1.data_;
	data_dstz = dstz.data_;
	data_dstz1 = dstz1.data_;
	data_srcx1 = srcx1.data_;
	data_srcy1 = srcy1.data_;
	data_srcz1 = srcz1.data_;
	data_src2 = src2.data_;

	int j, k;
	int stride_x2 = dstx.stride_[2];
	int stride_x1 = dstx.stride_[1];
	int stride_x0 = dstx.stride_[0];
	int stride_y2 = dsty.stride_[2];
	int stride_y1 = dsty.stride_[1];
	int stride_y0 = dsty.stride_[0];
	int stride_z2 = dstz.stride_[2];
	int stride_z1 = dstz.stride_[1];
	int stride_z0 = dstz.stride_[0];

	k = dstx.storage_.base(2);
	int kend = k + dstx.length_[2];
	for (; k < kend; k++) {
		RFloat *base_dstx = data_dstx + k * stride_x2;
		RFloat *base_dstx1 = data_dstx1 + k * stride_x2;
		RFloat *base_dsty = data_dsty + k * stride_y2;
		RFloat *base_dsty1 = data_dsty1 + k * stride_y2;
		RFloat *base_dstz = data_dstz + k * stride_z2;
		RFloat *base_dstz1 = data_dstz1 + k * stride_z2;
		RFloat *base_srcx1 = data_srcx1 + k * stride_x2;
		RFloat *base_srcx2 = data_src2 + k * stride_x2;
		RFloat *base_srcy1 = data_srcy1 + k * stride_y2;
		RFloat *base_srcy2 = data_src2 + k * stride_y2;
		RFloat *base_srcz1 = data_srcz1 + k * stride_z2;
		RFloat *base_srcz2 = data_src2 + k * stride_z2;

		j = dstx.storage_.base(1);
		int jend = j + dstx.length_[1];
		for (; j < jend; j++) {
			RFloat * base_dstx_1 = base_dstx + j * stride_x1;
			RFloat * base_dstx1_1 = base_dstx1 + j * stride_x1;
			RFloat * base_dsty_1 = base_dsty + j * stride_y1;
			RFloat * base_dsty1_1 = base_dsty1 + j * stride_y1;
			RFloat * base_dstz_1 = base_dstz + j * stride_z1;
			RFloat * base_dstz1_1 = base_dstz1 + j * stride_z1;
			RFloat * base_srcx1_1 = base_srcx1 + j * stride_x1;
			RFloat * base_srcx2_1 = base_srcx2 + j * stride_x1;
			RFloat * base_srcy1_1 = base_srcy1 + j * stride_y1;
			RFloat * base_srcz1_1 = base_srcz1 + j * stride_z1;

			int ist = dstx.storage_.base(0);
			int iend = ist + dstx.length_[0];

			local_update_msubadd_xyz(ist, iend, base_dstx_1, base_dstx1_1, base_srcx1_1, base_srcx2_1,base_dsty_1, base_dsty1_1, base_srcy1_1,base_dstz_1, base_dstz1_1, base_srcz1_1);
		}
	}


}



void local_DqdxAddWorksxyzMultiQ_finit(RDouble * __restrict dqdx, RDouble * __restrict dqdy, RDouble * __restrict dqdz, RDouble * __restrict worksx, RDouble * __restrict worksy, RDouble * __restrict worksz, RDouble * __restrict q, int ist, int ied) {
	for (int i = ist; i <= ied; i++) {
		RDouble qi = q[i];
		dqdx[i] = -worksx[i] * qi;
		dqdy[i] = -worksy[i] * qi;
		dqdz[i] = -worksz[i] * qi;
	}
}

void local_DqdxAddWorksxyzMultiQ_assign(int ifirst, int ilast,
		/*assign-*/RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, 
		/*assign+ or += */RDouble *dqdxl, RDouble *dqdyl, RDouble *dqdzl,
		const RDouble *worksx, const RDouble *worksy, const RDouble *worksz, const RDouble *q) {

	int i = ifirst;
	for (; i <= ilast; i++) {
		RDouble tmpx = worksx[i] * q[i];
		RDouble tmpy = worksy[i] * q[i];
		RDouble tmpz = worksz[i] * q[i];

		dqdx[i] = -tmpx; dqdxl[i] = tmpx; 
		dqdy[i] = -tmpy; dqdyl[i] = tmpy; 
		dqdz[i] = -tmpz; dqdzl[i] = tmpz; 
	}
}

void local_DqdxAddWorksxyzMultiQ_hybrid(int ifirst, int ilast,
		/*assign-*/RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, 
		/*assign+ or += */RDouble *dqdxl, RDouble *dqdyl, RDouble *dqdzl,
		const RDouble *worksx, const RDouble *worksy, const RDouble *worksz, const RDouble *q) {

	int i = ifirst;

	#pragma ivdep
	#pragma vector always
	for (; i <= ilast; i++) {
		RDouble tmp = q[i];
		RDouble tmpx = worksx[i] * tmp;
		RDouble tmpy = worksy[i] * tmp;
		RDouble tmpz = worksz[i] * tmp;

		dqdx[i] = -tmpx; dqdxl[i] += tmpx; 
		dqdy[i] = -tmpy; dqdyl[i] += tmpy; 
		dqdz[i] = -tmpz; dqdzl[i] += tmpz; 
	}
}

void local_update_Vanadd(RDouble MUSCLCoefXb, int ist, int iend, RDouble * base_dst_1, RDouble * base_dst1_1, RDouble * __restrict base_src1_1, RDouble * __restrict base_src2_1) {
            	int i;
 		#pragma ivdep
		#pragma vector always
		for (i = ist; i < iend; i++) {
			RDouble b2 = base_src2_1[i];
			RDouble b1 = base_src1_1[i];
			RDouble by = MUSCLCoefXb * b2;
			RDouble bx = MUSCLCoefXb * b1;
			const RDouble eps=1.0e-15;
			RDouble s1=1;
			RDouble s2=1;
			RDouble s1_1=1;
			RDouble s2_1=1;
			if(b1<0) s1=-1;
			if(by<0) s2=-1;
			if(b2<0) s1_1=-1;
			if(bx<0) s2_1=-1;
			base_dst_1[i]=(s1+s2)*b1*by/(ABS(b1+by)+eps);
			base_dst1_1[i]=(s1_1+s2_1)*b2*bx/(ABS(b2+bx)+eps);
			}

}

void MinModOpt(RDouble value,int ifirst,int iend,RDouble *dst1,RDouble *dst2,RDouble *__restrict src1,RDouble *__restrict srcsl1,RDouble * __restrict srcpl1)
{
	int i;
#pragma ivdep
#pragma vector always
	for(i = ifirst; i < iend;i++)
	{
		RDouble src2 = src1[i] - srcsl1[i];
		RDouble src3 = srcpl1[i] - src1[i];
		RDouble s1=1;
		RDouble s2=1;
		RDouble min;

		if(src2<0) s1=-1;
		if(src3<0) s2=-1;
		
		min = (ABS(src2) < ABS(src3)) ? ABS(src2):ABS(src3);
		RDouble s3 = (0.5*( s1 + s2 ) * min)*value;

		dst1[i]=src1[i] + s3;
		dst2[i]=src1[i] - s3;
	}

}
void VANALBADAOpt(RDouble value,int ifirst,int iend,RDouble *dst1,RDouble *dst2,RDouble *__restrict src1,RDouble *__restrict srcsl1,RDouble * __restrict srcpl1)
{
    const RDouble eps = 1.0e-7;
	RDouble src2;
	RDouble src3;
    RFloat cx,cy,cd;
    int i;

   #pragma ivdep
   #pragma vector always
	for(i = ifirst; i < iend;i++)
	{
		src2 = src1[i] - srcsl1[i];
		src3 = srcpl1[i] - src1[i];
		cx = src2 * src2 + eps;
		cy = src3 * src3 + eps;
	    cd = 1.0 / ( cx + cy );

	    dst1[i]=src1[i]+value*(0.5 * ( SIGN(1.0, src2) + SIGN(1.0, src3) ) * ( src2 * src3 * ABS( src2 + src3 ) ) * cd); 
		dst2[i]=src1[i]-value*(0.5 * ( SIGN(1.0, src2) + SIGN(1.0, src3) ) * ( src2 * src3 * ABS( src2 + src3 ) ) * cd); 
	}
}

void VanalbadaCLZOpt(RDouble value,int ifirst,int iend,RDouble *dst1,RDouble *dst2,RDouble *__restrict src1,RDouble *__restrict srcsl1,RDouble * __restrict srcpl1)
{
    const RDouble eps = 1.0e-7;
	RDouble src2;
	RDouble src3;
    RFloat cx,cy,cd;
    int i;
   #pragma ivdep
   #pragma vector always
	for(i = ifirst; i < iend;i++)
	{
		src2 = src1[i] - srcsl1[i];
		src3 = srcpl1[i] - src1[i];
		cx = src2 * src2 + eps;
		cy = src3 * src3 + eps;
	    cd = 1.0 / ( cx + cy );

		dst1[i]=src1[i]+value*(0.5 * ( SIGN(1.0, src2) + SIGN(1.0, src3) ) * ( cx * src3 + cy * src2 ) * cd); 
		dst2[i]=src1[i]-value*(0.5 * ( SIGN(1.0, src2) + SIGN(1.0, src3) ) * ( cx * src3 + cy * src2 ) * cd); 
	}
}
void VanleerOpt(RDouble value,int ifirst,int iend,RDouble *dst1,RDouble *dst2,RDouble *__restrict src1,RDouble *__restrict srcsl1,RDouble * __restrict srcpl1)
{
    const RDouble eps = 1.0e-15;
	RDouble src2;
	RDouble src3; 
    int i;
   #pragma ivdep
   #pragma vector always
	for(i = ifirst; i < iend;i++)
	{
		src2 = src1[i] - srcsl1[i];
		src3 = srcpl1[i] - src1[i];

		dst1[i]=src1[i]+value*(( SIGN(1.0, src2) + SIGN(1.0, src3) ) * src2 * src3 / ( ABS( src2 + src3 ) + eps )); 
		dst2[i]=src1[i]-value*(( SIGN(1.0, src2) + SIGN(1.0, src3) ) * src2 * src3 / ( ABS( src2 + src3 ) + eps )); 
	}
}
void SmoothOpt(RDouble value,int ifirst,int iend,RDouble *dst1,RDouble *dst2,RDouble *__restrict src1,RDouble *__restrict srcsl1,RDouble * __restrict srcpl1)
{
    const RDouble eps = 1.0e-15;
	RDouble src2;
	RDouble src3; 
    int i;

	RFloat Smooth1;
    RFloat Smooth2;

   #pragma ivdep
   #pragma vector always
	for(i = ifirst; i < iend;i++)
	{
		src2 = src1[i] - srcsl1[i];
		src3 = srcpl1[i] - src1[i];
	
		Smooth1= ( 2. * src2 * src3 + eps ) / (  src2 * src2 + src3 * src3 + eps );
	    Smooth2= ((1.+Smooth1)*src2+ (1.-Smooth1)*src3)*Smooth1/2.;   //0.5 * ( SIGN(1.0, x) + SIGN(1.0, y) )*

		dst1[i]=src1[i]+value*Smooth2; 
		dst2[i]=src1[i]-value*Smooth2; 
	}
}

void MinvanOpt(RDouble value,int ifirst,int iend,RDouble *dst1,RDouble *dst2,RDouble *__restrict src1,RDouble *__restrict srcsl1,RDouble * __restrict srcpl1)
{
    const RDouble eps = 1.0e-15;
	RDouble src2;
	RDouble src3; 
    int i;

	RFloat z,z2,z3,xy2,absx,absy,min_van;

   #pragma ivdep
   #pragma vector always
	for(i = ifirst; i < iend;i++)
	{
		src2 = src1[i] - srcsl1[i];
		src3 = srcpl1[i] - src1[i];

		z = 0.5 * (SIGN(1.0, src2) + SIGN(1.0, src3) );
		absx = abs(src2);
		absy = abs(src3);
		xy2 = (src2+src2)*src3;
		z2 =  xy2/( src2*src2 + src3*src3 + eps);
		z3 =  xy2/( absx + absy + eps );
		min_van = z * ( (1.0-z2)* min( absx , absy ) + z2 * z3 );

		dst1[i]=src1[i]+value*min_van; 
		dst2[i]=src1[i]-value*min_van; 
	}
}
void Q_S_Q_Opt(RDouble value,int ifirst,int iend,RDouble *dst1,RDouble *dst2,RDouble *__restrict src1,RDouble *__restrict srcsl1,RDouble * __restrict srcpl1)
{
	RDouble src2;
	RDouble src3; 
    int i;

   #pragma ivdep
   #pragma vector always
	for(i = ifirst; i < iend;i++)
	{
		src2 = src1[i] - srcsl1[i];
		src3 = srcpl1[i] - src1[i];

		dst1[i]=src1[i]+value*src2; 
		dst2[i]=src1[i]-value*src3; 
	}
}

void ResComputing(PHArray<RFloat, 4> dst,const PHArray<RFloat, 4> &src,int il1,int jl1,int kl1)
{
	#pragma ivdep
	#pragma vector always	

	int i,j, k, m,iend,jend,kend,mend;

    declare_data(dst); 
    declare_data(src); 
 	 	
	declare_stride_4d(src);

	int J_p_jl1;
    int K_p_kl1;

   m = src.storage_.base(3);
   mend = m + src.length_[3]; 
   for (; m< mend; m++)
    {   
		declare_basex(dst,3,data,m);
		declare_base_suffixx(src,3,data,m,m,);		
		k = src.storage_.base(2);
		kend = k + src.length_[2];
		for (  ; k<kend; k++)
		{
			K_p_kl1=k+kl1;
			declare_basex(dst,2,3,k);
			declare_base_suffixx(src,2,3,K_p_kl1,k_p_src,m);	
			declare_base_suffixx(src,2,3,k,src_k,m);
			j = src.storage_.base(1);
			jend = j + src.length_[1];	
			for ( ;j< jend; j++)
			{
				J_p_jl1=j+jl1;
				declare_basex(dst,1,2,j);
				declare_base_suffixx(src,1,2,J_p_jl1,j_p_src,k_p_src);
				declare_base_suffixx(src,1,2,j,src_j,src_k);               
				 i = src.storage_.base(0);
			     iend = i + src.length_[0];	
				 for (; i<iend ; i++)
				 {	
					declare_base_namex(dst,1)[i]-=(declare_base_name_suffixx(src,1,j_p_src)+il1)[i]-declare_base_name_suffixx(src,1,src_j)[i];	

				 } 
				  
			}
		}
	}
}

void AssignInPHOpt(PHArray<RFloat, 3> dst1,PHArray<RFloat, 3> dst2,PHArray<RFloat, 3> dst3,PHArray<RFloat, 3> dst4,PHArray<RFloat, 3> dst5,PHArray<RFloat, 3> dst6,PHArray<RFloat, 3> dst7,PHArray<RFloat, 3> dst8,PHArray<RFloat, 3> dst9,const PHArray<RFloat, 3> &src1,const PHArray<RFloat, 3> &src2,const PHArray<RFloat, 3> &src3,const PHArray<RFloat, 3> &src4,const PHArray<RFloat, 3> &src5,const PHArray<RFloat, 3> &src6,const PHArray<RFloat, 3> &src7,const PHArray<RFloat, 3> &src8,const PHArray<RFloat, 3> &src9)
{
	 declare_data(dst1); 
     declare_data(dst2); 
	 declare_data(dst3);
	 declare_data(dst4);
	 declare_data(dst5); 
	 declare_data(dst6);
	 declare_data(dst7);
	 declare_data(dst8);
	 declare_data(dst9); 
	 declare_data(src1); 
     declare_data(src2); 
	 declare_data(src3);
	 declare_data(src4);
	 declare_data(src5); 
	 declare_data(src6);
	 declare_data(src7);
	 declare_data(src8);
	 declare_data(src9); 
 	
	 declare_stride_3d(src1);

    int i,j, k, iend,jend,kend;
	k = src1.storage_.base(2);
	kend = k + src1.length_[2];
	for (  ; k<kend; k++)
	{			 
		declare_basex(dst1,2,data,k);
		declare_basex(dst2,2,data,k);
		declare_basex(dst3,2,data,k);
		declare_basex(dst4,2,data,k);
		declare_basex(dst5,2,data,k);
		declare_basex(dst6,2,data,k);
		declare_basex(dst7,2,data,k);
		declare_basex(dst8,2,data,k);
		declare_basex(dst9,2,data,k);
		declare_basex(src1,2,data,k);
		declare_basex(src2,2,data,k);
		declare_basex(src3,2,data,k);
		declare_basex(src4,2,data,k);
		declare_basex(src5,2,data,k);
		declare_basex(src6,2,data,k);
		declare_basex(src7,2,data,k);
		declare_basex(src8,2,data,k);
		declare_basex(src9,2,data,k);
		j = src1.storage_.base(1);
		jend = j + src1.length_[1];	
		for ( ;j< jend; j++)
		{
            declare_basex(dst1,1,2,j);
			declare_basex(dst2,1,2,j);
			declare_basex(dst3,1,2,j);
			declare_basex(dst4,1,2,j);
			declare_basex(dst5,1,2,j);
			declare_basex(dst6,1,2,j);
			declare_basex(dst7,1,2,j);
			declare_basex(dst8,1,2,j);
			declare_basex(dst9,1,2,j); 
		    declare_basex(src1,1,2,j);
			declare_basex(src2,1,2,j);
			declare_basex(src3,1,2,j);
			declare_basex(src4,1,2,j);
			declare_basex(src5,1,2,j);
			declare_basex(src6,1,2,j);
			declare_basex(src7,1,2,j);
			declare_basex(src8,1,2,j);
			declare_basex(src9,1,2,j);  
			i = src1.storage_.base(0);
			iend = i + src1.length_[0];	
			for (; i<iend ; i++)
			{		 		
				declare_base_namex(dst1,1)[i] = declare_base_namex(src1,1)[i];
				declare_base_namex(dst2,1)[i] = declare_base_namex(src2,1)[i];
				declare_base_namex(dst3,1)[i] = declare_base_namex(src3,1)[i];

				declare_base_namex(dst4,1)[i] = declare_base_namex(src4,1)[i];
				declare_base_namex(dst5,1)[i] = declare_base_namex(src5,1)[i];
				declare_base_namex(dst6,1)[i] = declare_base_namex(src6,1)[i];

				declare_base_namex(dst7,1)[i] = declare_base_namex(src7,1)[i];
				declare_base_namex(dst8,1)[i] = declare_base_namex(src8,1)[i];
				declare_base_namex(dst9,1)[i] = declare_base_namex(src9,1)[i];
			} 
				  
		}
	}
	

}

void ZeroAssig(PHArray<RFloat, 4> &src)
{
	int i,j, k, m,mend,kend,jend,iend;
	declare_data(src);
	declare_stride_4d(src);
	m = src.storage_.base(3);
    mend = m + src.length_[3];

	for (; m < mend; m++)
	{
		declare_basex(src,3,data,m);
		k = src.storage_.base(2);
	    kend = k + src.length_[2];
		for (; k<kend;k++)
		{
			declare_basex(src,2,3,k);	
			j = src.storage_.base(1);
			jend = j + src.length_[1];
			for (; j<jend;j++)
			{
			   declare_basex(src,1,2,j);
			   	i = src.storage_.base(0);
				iend = i + src.length_[0];	
			   for (; i<iend; i++)
			   {
				   declare_base_namex(src,1)[i]=0.0;
			   }
			}
		}
	}
}

void ParaVanleer(PHArray<RFloat, 4> dst, PHArray<RFloat, 4> dst1, const PHArray<RFloat, 4> &src1, const PHArray<RFloat, 4> &src2, RDouble MUSCLCoefXb) {
	RFloat * data_dst, *data_dst1, *data_src1, *data_src2;

	data_dst = dst.data_;
	data_dst1 = dst1.data_;
	data_src1 = src1.data_;
	data_src2 = src2.data_;

	int j, k, m;
	int stride3 = dst.stride_[3];
	int stride2 = dst.stride_[2];
	int stride1 = dst.stride_[1];
	int stride0 = dst.stride_[0];
	int stride31 = dst1.stride_[3];
	int stride21 = dst1.stride_[2];
	int stride11 = dst1.stride_[1];
	int stride01 = dst1.stride_[0];
	
	m = dst.storage_.base(3);
	int mend = m + dst.length_[3];
	for (; m < mend; m++) {
		RFloat *base_dst = data_dst + m * stride3;
		RFloat *base_dst1 = data_dst1 + m * stride3;
		RFloat *base_src1 = data_src1 + m * stride3;
		RFloat *base_src2 = data_src2 + m * stride3;

	k = dst.storage_.base(2);
	int kend = k + dst.length_[2];
	for (; k < kend; k++) {
		RFloat *base_dst_0 = base_dst + k * stride2;
		RFloat *base_dst1_0 = base_dst1 + k * stride2;
		RFloat *base_src1_0 = base_src1 + k * stride2;
		RFloat *base_src2_0 = base_src2 + k * stride2;

		j = dst.storage_.base(1);
		int jend = j + dst.length_[1];
		for (; j < jend; j++) {
			RFloat * base_dst_1 = base_dst_0 + j * stride1;
			RFloat * base_dst1_1 = base_dst1_0 + j * stride1;
			RFloat * base_src1_1 = base_src1_0 + j * stride1;
			RFloat * base_src2_1 = base_src2_0 + j * stride1;

			int ist = dst.storage_.base(0);
			int iend = ist + dst.length_[0];

			local_update_Vanadd(MUSCLCoefXb, ist, iend, base_dst_1, base_dst1_1, base_src1_1, base_src2_1);

		}
	}
	}
}
}
#endif

