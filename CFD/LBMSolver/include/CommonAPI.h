
#pragma once
#include "DevEnv.h"
#pragma warning(disable:4251)

#ifndef COMMONAPI	// 以项目名命名的接口控制宏
#ifdef _WIN32
#ifdef COMMONLIBDLL  
#define COMMONAPI _declspec(dllexport)  
#else
#define COMMONAPI _declspec(dllimport)  
#endif 
#else
#define COMMONAPI // linux系统不需要
#endif
#endif
//
//#define SOLUTION_DIR "D:\\workspace\\nnwlbm\\lbm-dev\\"
//#define JOINT(p,f)   p##f 
//
//#ifndef COMMONLIBDLL
//#ifdef _WIN32
//#ifdef _WIN64
//#ifdef _DEBUG
//#pragma comment(lib, JOINT(SOLUTION_DIR, "x64\\Debug\\common.lib"))
//#else //_DEBUG
//#pragma comment(lib, JOINT(SOLUTION_DIR, "x64\\Release\\common.lib"))
//#endif // _DEBUG
//#else // _WIN64
//#ifdef _DEBUG
//#pragma comment(lib, JOINT(SOLUTION_DIR, "Debug\\common.lib"))
//#else //_DEBUG
//#pragma comment(lib, JOINT(SOLUTION_DIR, "Release\\common.lib"))
//#endif // _DEBUG
//#endif // _WIN64
//#endif // _WIN32
//#endif // _COMMONLIBDLL
