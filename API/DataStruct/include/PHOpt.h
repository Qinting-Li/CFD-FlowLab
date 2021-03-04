//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      PHOpt.h
//! @brief     Define some optimization tools defined by PARATERA.
//! @author    PARATERA.

#pragma once
namespace PHSPACE
{
#ifdef USE_ARRAY_OPT
#define declare_stride_1d(fyarray) int stride_0 = fyarray.get_stride(0) //!declare stride_0 = q.stride_[1];stride_1 = q.stride_[2];stride_2 = q.stride_[3];stride_3 = q.stride_[4];
#define declare_stride_2d(fyarray) int stride_1 = fyarray.get_stride(1); declare_stride_1d(fyarray)
#define declare_stride_3d(fyarray) int stride_2 = fyarray.get_stride(2); declare_stride_2d(fyarray)
#define declare_stride_4d(fyarray) int stride_3 = fyarray.get_stride(3); declare_stride_3d(fyarray)
#define declare_stride_5d(fyarray) int stride_4 = fyarray.get_stride(4); declare_stride_4d(fyarray)
#define declare_data(fyarray) RFloat *fyarray##_base##_data = fyarray.getData()//!eg.*gama_base_data = gama.data_
#define declare_basex(fyarray, cur, pre, indx) RFloat *fyarray##_base_##cur = fyarray##_base_##pre + (indx)*stride_##cur//!solve the length of fourth dimention; declare *fyarray_base_cur = fyarray_base_pre + indx * stride_cur;
#define declare_base_namex(fyarray, cur) fyarray##_base_##cur//!fyarray_base_cur[i];
#define declare_base_suffixx(fyarray, cur, pre, indx, cur_suf, pre_suf) RFloat *fyarray##cur_suf##_base_##cur = fyarray##pre_suf##_base_##pre + (indx)*stride_##cur//!solve the length of fourth dimention; declare *fyarray_cur_suf_base_cur = fyarray_pre_suf_base_pre + indx * stride_cur;
#define declare_base_name_suffixx(fyarray, cur, suffix) fyarray##suffix##_base_##cur//!fyarray_suffix_base_cur[i];
#endif
}