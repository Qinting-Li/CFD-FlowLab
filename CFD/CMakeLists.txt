#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +  
#          P   P  H   H  E      NN   N  G      L      E       I          +
#          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
#          P      H   H  E      N  N N  G   G  L      E       I          +
#          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
#------------------------------------------------------------------------+
#          Platform for Hybrid Engineering Simulation of Flows           +
#           China Aerodynamics Research and Development Center           +
#                     (C) Copyright, Since 2010                          +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! @file      CMakeLists.txt
#! @brief     It defines the project.
#! @author    Bell, Xu Qingxin, He Kun, He Xin.
cmake_minimum_required(VERSION 2.6)

include(PHMacro.cmake)

#set working directory:
set(EXEC_WORKING_DIR PHengLEIv3d0)

#set cmake environment:
set(CMAKE_USE_RELATIVE_PATHS ON)

#output message:
message(STATUS "Welcome to Build Project ${CMAKE_PROJECT_NAME}.")
message(STATUS "The OS is ${CMAKE_SYSTEM_NAME}, version ${CMAKE_SYSTEM_VERSION}.")
message(STATUS "The working directory is ${CMAKE_BINARY_DIR}/${EXEC_WORKING_DIR}/${DLL_WORKING_DIR}.")
message(STATUS "--- Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

# Define some pre-definition.
include(PHMacro.cmake)
if(USE_OCC)
  message(STATUS "USE_OCC is turned on!")
  add_definitions(-DUSE_OCC)
endif()
if(USE_DOUBLE_PRECISION)
  message(STATUS "USE_DOUBLE_PRECISION is turned on!")
  add_definitions(-DPH_DOUBLE_PREC)
endif()
if(USE_SMARTARRAY)
  message(STATUS "USE_SMARTARRAY is turned on!")
  add_definitions(-DUSE_SMARTARRAY)
endif()
if(USE_ARRAY_OPT)
  message(STATUS "USE_ARRAY_OPT is turned on!")
  add_definitions(-DUSE_ARRAY_OPT)
endif()
if(WIN32)
if(USE_Windows_X64)
  message(STATUS "USE_Windows_X64 is turned on!")
  add_definitions(-DUSE_WINDOWS_X64)
endif(USE_Windows_X64)
endif(WIN32)
#1. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory <CFD>, Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

set(CFD_TARGET_NAME "PHCFD")
set(LIBRARY_OUTPUT_PATH)
set(LIB_DataStruct_EXPORT_PATH)
set(LIB_Data_EXPORT_PATH)
set(LIB_Math_EXPORT_PATH)
set(LIB_PHMPI_EXPORT_PATH)
set(LIB_Toolkit_EXPORT_PATH)
set(LIB_IO_EXPORT_PATH)
set(LIB_Geometry_EXPORT_PATH)
set(LIB_Physics_EXPORT_PATH)
set(LIB_PreProcess_EXPORT_PATH)
set(LIB_PostProcess_EXPORT_PATH)
set(LIB_Common_EXPORT_PATH)
if(WIN32)

if(USE_Windows_X64)
 set(LIB_DataStruct_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/DataStruct/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIB_Data_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Data/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIB_Math_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Math/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIB_PHMPI_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/FYMPI/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIB_Toolkit_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Toolkit/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIB_IO_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/IO/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIB_Geometry_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Geometry/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIB_Physics_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Physics/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIB_PreProcess_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/PreProcess/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIB_PostProcess_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/PostProcess/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIB_Common_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Common/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(CMAKE_EXE_LINKER_FLAGS "/machine:x64")
else(USE_Windows_X64)
  set(LIB_DataStruct_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/DataStruct/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_Data_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Data/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_Math_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Math/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_PHMPI_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/FYMPI/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_Toolkit_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Toolkit/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_IO_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/IO/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_Geometry_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Geometry/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_Physics_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Physics/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
 set(LIB_PreProcess_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/PreProcess/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
 set(LIB_PostProcess_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/PostProcess/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
 set(LIB_Common_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Common/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(USE_Windows_X64)
else(WIN32)
  set(LIB_DataStruct_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/DataStruct/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_Data_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Data/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_Math_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Math/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_PHMPI_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/FYMPI/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_Toolkit_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Toolkit/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_IO_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/IO/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_Geometry_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Geometry/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_Physics_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Physics/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_PreProcess_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/PreProcess/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_PostProcess_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/PostProcess/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_Common_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/Common/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)
link_directories("${LIB_DataStruct_EXPORT_PATH}")
link_directories("${LIB_Data_EXPORT_PATH}")
link_directories("${LIB_Math_EXPORT_PATH}")
link_directories("${LIB_PHMPI_EXPORT_PATH}")
link_directories("${LIB_Toolkit_EXPORT_PATH}")
link_directories("${LIB_IO_EXPORT_PATH}")
link_directories("${LIB_Geometry_EXPORT_PATH}")
link_directories("${LIB_Physics_EXPORT_PATH}")
link_directories("${LIB_PreProcess_EXPORT_PATH}")
link_directories("${LIB_PostProcess_EXPORT_PATH}")
link_directories("${LIB_Common_EXPORT_PATH}")
message(STATUS "CFD library path is: ${LIBRARY_OUTPUT_PATH}.")

#set options
if (NOT BUILD_LIBRARY_TYPE)
  set (BUILD_LIBRARY_TYPE "Static" CACHE STRING "The type of PHengLEI libraries" FORCE)
  SET_PROPERTY(CACHE BUILD_LIBRARY_TYPE PROPERTY STRINGS Shared Static)
endif()

#3. set environment for MPICH and 3rdParty lib directory:
if(WIN32)
  set(MSMPI_HOME_DIRECTORY "C:/Program Files (x86)/Microsoft SDKs/MPI" CACHE PATH "path to installed MPICH")
  include_directories("${MSMPI_HOME_DIRECTORY}/Include")
  if(USE_Windows_X64)
    link_directories("${MSMPI_HOME_DIRECTORY}/Lib/x64")
    set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x64")
  else(USE_Windows_X64) 
    link_directories("${MSMPI_HOME_DIRECTORY}/Lib/x86")
    set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x86")
  endif(USE_Windows_X64)
else(WIN32)
  set(CMAKE_CXX_COMPILER  mpicxx)
  set(CMAKE_C_COMPILER    mpicc)
  
  set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/3rdparty/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif(WIN32)
# specify directory in which the linker will look for libraries for myself libraries:
message(STATUS "The 3rdparty path is: ${THIRDPARTY_LIB_DIR}.") 
link_directories("${THIRDPARTY_LIB_DIR}")

set(DataStruct_H_PATH   ${PROJECT_SOURCE_DIR}/../API/DataStruct/include)
set(Data_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Data/include)
set(Math_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Math/include)
set(PHMPI_H_PATH   ${PROJECT_SOURCE_DIR}/../API/FYMPI/include)
set(Toolkit_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Toolkit/include)
set(IO_H_PATH   ${PROJECT_SOURCE_DIR}/../API/IO/include)
set(Geometry_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Geometry/include)
set(Physics_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Physics/include)
set(PreProcess_H_PATH   ${PROJECT_SOURCE_DIR}/../API/PreProcess/include)
set(PostProcess_H_PATH   ${PROJECT_SOURCE_DIR}/../API/PostProcess/include)
set(Common_H_PATH   ${PROJECT_SOURCE_DIR}/../API/Common/include)
include_directories("${DataStruct_H_PATH}")
include_directories("${Data_H_PATH}")
include_directories("${Math_H_PATH}")
include_directories("${PHMPI_H_PATH}")
include_directories("${Toolkit_H_PATH}")
include_directories("${IO_H_PATH}")
include_directories("${Geometry_H_PATH}")
include_directories("${Physics_H_PATH}")
include_directories("${PreProcess_H_PATH}")
include_directories("${PostProcess_H_PATH}")
include_directories("${Common_H_PATH}")
set(HDF5_H_PATH    ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/include)
include_directories("${HDF5_H_PATH}")
set(TECIO_H_PATH   ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/TECIO/include)
include_directories("${TECIO_H_PATH}")
set(MEITIS_H_PATH  ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/metis/include)
include_directories("${MEITIS_H_PATH}")
set(CGNS_H_PATH    ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/cgns/include)
include_directories("${CGNS_H_PATH}")
#set(API_H_DIR ${PROJECT_SOURCE_DIR}/3rdparty)

#5. set current directory:
# set directory of Solver
set(SOLVER_H_PATH    Solver/include)
set(SOLVER_H_LIST
${SOLVER_H_PATH}/CFDSolver.h
${SOLVER_H_PATH}/Controller.h
${SOLVER_H_PATH}/Param_CFDSolver.h
)
source_group(Solver\\include FILES ${SOLVER_H_LIST}) 
include_directories(${SOLVER_H_PATH})

set(SOLVER_CPP_PATH      Solver/src)
set(SOLVER_CPP_LIST
${SOLVER_CPP_PATH}/CFDSolver.cpp
${SOLVER_CPP_PATH}/Controller.cpp
${SOLVER_CPP_PATH}/Param_CFDSolver.cpp
)
source_group(Solver\\src FILES ${SOLVER_CPP_LIST})

# set directory of NSSolver:
set(NSSOLVER_H_PATH    NSSolver/include)
set(NSSOLVER_H_LIST
${NSSOLVER_H_PATH}/NSSolver.h
${NSSOLVER_H_PATH}/NSSolverStruct.h
${NSSOLVER_H_PATH}/NSSolverStructFD.h
${NSSOLVER_H_PATH}/NSSolverStructFV.h
${NSSOLVER_H_PATH}/NSSolverUnstruct.h
${NSSOLVER_H_PATH}/Boundary.h
${NSSOLVER_H_PATH}/MultiGridOperation.h
${NSSOLVER_H_PATH}/Param_NSSolver.h
${NSSOLVER_H_PATH}/Param_NSSolverUnstruct.h
${NSSOLVER_H_PATH}/Param_NSSolverStruct.h
)
source_group(NSSolver\\include FILES ${NSSOLVER_H_LIST}) 
include_directories(${NSSOLVER_H_PATH})

set(NSSOLVER_CPP_PATH      NSSolver/src)
set(NSSOLVER_CPP_LIST
${NSSOLVER_CPP_PATH}/NSSolver.cpp
${NSSOLVER_CPP_PATH}/NSSolverStruct.cpp
${NSSOLVER_CPP_PATH}/NSSolverStructFD.cpp
${NSSOLVER_CPP_PATH}/NSSolverStructFV.cpp
${NSSOLVER_CPP_PATH}/NSSolverUnstruct.cpp
${NSSOLVER_CPP_PATH}/Boundary.cpp
${NSSOLVER_CPP_PATH}/Param_NSSolver.cpp
${NSSOLVER_CPP_PATH}/Param_NSSolverUnstruct.cpp
${NSSOLVER_CPP_PATH}/Param_NSSolverStruct.cpp
)
source_group(NSSolver\\src FILES ${NSSOLVER_CPP_LIST})

# set directory of TurbSolver:
set(TurbSOLVER_H_PATH        TurbSolver/include)
set(TurbSOLVER_H_LIST
${TurbSOLVER_H_PATH}/TurbSolver.h
${TurbSOLVER_H_PATH}/TurbSolverStruct.h
${TurbSOLVER_H_PATH}/TurbSolverUnstr.h
${TurbSOLVER_H_PATH}/Param_TurbSolver.h
${TurbSOLVER_H_PATH}/Param_TurbSolverStruct.h
${TurbSOLVER_H_PATH}/Param_TurbSolverUnstruct.h
)
source_group(TurbSolver\\include FILES ${TurbSOLVER_H_LIST}) 
include_directories(${TurbSOLVER_H_PATH})

set(TurbSOLVER_CPP_PATH      TurbSolver/src)
set(TurbSOLVER_CPP_LIST
${TurbSOLVER_CPP_PATH}/TurbSolver.cpp
${TurbSOLVER_CPP_PATH}/TurbSolverStruct.cpp
${TurbSOLVER_CPP_PATH}/TurbSolverUnstr.cpp
${TurbSOLVER_CPP_PATH}/Param_TurbSolver.cpp
${TurbSOLVER_CPP_PATH}/Param_TurbSolverStruct.cpp
${TurbSOLVER_CPP_PATH}/Param_TurbSolverUnstruct.cpp
)
source_group(TurbSolver\\src FILES ${TurbSOLVER_CPP_LIST})

# set directory of HOSolverUnstruct:
set(HOSOLVERUNSTRUCT_H_PATH        HOSolverUnstruct/include)
set(HOSOLVERUNSTRUCT_H_LIST
${HOSOLVERUNSTRUCT_H_PATH}/HOBasisFunction.h
${HOSOLVERUNSTRUCT_H_PATH}/HOBoundary.h
${HOSOLVERUNSTRUCT_H_PATH}/HODefine.h
${HOSOLVERUNSTRUCT_H_PATH}/HOGaussJacobiQuadrature.h
${HOSOLVERUNSTRUCT_H_PATH}/HOGeometryStructure.h
${HOSOLVERUNSTRUCT_H_PATH}/HOInviscidFlux.h
${HOSOLVERUNSTRUCT_H_PATH}/HOShapeFunction.h
${HOSOLVERUNSTRUCT_H_PATH}/HOSimulation.h
${HOSOLVERUNSTRUCT_H_PATH}/HOSolverUnstruct.h
${HOSOLVERUNSTRUCT_H_PATH}/HOSolverUnstructParam.h
${HOSOLVERUNSTRUCT_H_PATH}/HOStandardElement.h
)
source_group(HOSolverUnstruct\\include FILES ${HOSOLVERUNSTRUCT_H_LIST}) 
include_directories(${HOSOLVERUNSTRUCT_H_PATH})

set(HOSOLVERUNSTRUCT_CPP_PATH      HOSolverUnstruct/src)
set(HOSOLVERUNSTRUCT_CPP_LIST
${HOSOLVERUNSTRUCT_CPP_PATH}/HOBasisFunction.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOBoundary.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOGaussJacobiQuadrature.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOGeometryDtructure.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOInviscidFlux.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOShapeFunction.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOSimulation.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOSolverUnstruct.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOSolverUnstructParam.cpp
${HOSOLVERUNSTRUCT_CPP_PATH}/HOStandardElement.cpp
)
source_group(HOSolverUnstruct\\src FILES ${HOSOLVERUNSTRUCT_CPP_LIST})



#20.1 set header file list
set(CFD_HEAD_LIST 
${SOLVER_H_LIST}
${NSSOLVER_H_LIST}
${TurbSOLVER_H_LIST}
${HOSOLVERUNSTRUCT_H_LIST}
)

#20.2 set source file list:
set(CFD_CPP_LIST 
${SOLVER_CPP_LIST}
${NSSOLVER_CPP_LIST}
${TurbSOLVER_CPP_LIST}
${HOSOLVERUNSTRUCT_CPP_LIST}
)

#21. adds flags to the compiler command line for sources in the current directory and below:
add_definitions(-DPH_PARALLEL -DPH_CGNS)

# For OCC.
set(OCC_LIBS)
if(USE_OCC)
  message(STATUS "Using OpenCasCade ...... ")
  set(OCC_LIBS_REQUIRED
      # subset of DataExchange
      TKSTEP TKSTEP209 TKSTEPAttr TKSTEPBase TKIGES TKXSBase
      # ModelingAlgorithms
      TKOffset TKFeat TKFillet TKBool TKShHealing TKMesh TKHLR TKBO TKPrim
      TKTopAlgo TKGeomAlgo
      # ModelingData
      TKBRep TKGeomBase TKG3d TKG2d
      # FoundationClasses
      # TKAdvTools -- not necessary? (and removed from OCC 6.8)
      TKMath TKernel)
  
  # Set the path of the OCC library.
  set(3RDPARTY_LIB_PATH ${PROJECT_SOURCE_DIR}/3rdparty/occ690-vc11-bit32-release)
  set(OCC_LIB_PATH ${3RDPARTY_LIB_PATH} CACHE PATH "path to OpenCasCade library.")
  message(STATUS "Path to OpenCasCade library: " ${OCC_LIB_PATH})
  
  # Check if all the library exist -- it's need to be checked. 
  foreach(OCC ${OCC_LIBS_REQUIRED})     
	  list(APPEND OCC_LIBS ${OCC})
      unset(OCC CACHE)
  endforeach(OCC)
  
  if(WIN32)           
    MESSAGE(STATUS "Using Win32 OpenCasCade.")	
    MESSAGE(STATUS "${OCC_LIB_PATH}/lib")
    MESSAGE(STATUS "${OCC_LIB_PATH}/bin")	
    link_directories("${OCC_LIB_PATH}/lib")
    link_directories("${OCC_LIB_PATH}/bin") 
	include_directories("${OCC_LIB_PATH}/inc")
  else(WIN32)
    MESSAGE(STATUS "Using NON-Win32 OpenCasCade.")
    MESSAGE(STATUS "${OCC_LIB_PATH}/lin64/gcc/lib")	
    MESSAGE(STATUS "${OCC_LIB_PATH}/inc")	         
    link_directories("${OCC_LIB_PATH}/lin64/gcc/lib") 
    include_directories("${OCC_LIB_PATH}/inc")        
  endif(WIN32)
endif()

# Some pre-definitions that occ need.
if(USE_OCC)
  #target_link_libraries(${EXECUTABLE_TARGET_NAME} ${OCC_LIBS})  
  add_definitions(-DWNT)
  add_definitions(-DNOTCL)
endif()

# build executable program or library to the project using the specified source files and header files:
if ("${BUILD_LIBRARY_TYPE}" STREQUAL "Shared")
  # Building Shared Library.   
  message(STATUS "Building Shared Library ${CFD_TARGET_NAME} ...")
  add_library(${CFD_TARGET_NAME} SHARED ${CFD_CPP_LIST} ${CFD_HEAD_LIST})
elseif ("${BUILD_LIBRARY_TYPE}" STREQUAL "Static")
  # Building Static Library.
  message(STATUS "Building Static Library ${CFD_TARGET_NAME} ...")
  add_library(${CFD_TARGET_NAME} STATIC ${CFD_CPP_LIST} ${CFD_HEAD_LIST}) 
endif() 


# link API library:
set(DataStruct_TARGET_NAME "PHDataStruct")
set(Data_TARGET_NAME "PHData")
set(Math_TARGET_NAME "PHMath")
set(PHMPI_TARGET_NAME "PHMPI")
set(Toolkit_TARGET_NAME "PHToolkit")
set(IO_TARGET_NAME "PHIO")
set(Geometry_TARGET_NAME "PHGeometry")
set(Physics_TARGET_NAME "PHPhysics")
set(PreProcess_TARGET_NAME "PHPreProcess")
set(PostProcess_TARGET_NAME "PHPostProcess")
set(Common_TARGET_NAME "PHCommon")
list(APPEND LINK_LIBRARY_LIST ${Physics_TARGET_NAME} ${PreProcess_TARGET_NAME} ${PostProcess_TARGET_NAME} ${Common_TARGET_NAME} ${Geometry_TARGET_NAME} ${IO_TARGET_NAME} ${Toolkit_TARGET_NAME} ${PHMPI_TARGET_NAME} ${DataStruct_TARGET_NAME} ${Data_TARGET_NAME} ${Math_TARGET_NAME} metis parmetis cgns hdf5)

if(WIN32)
  list(APPEND LINK_LIBRARY_LIST msmpi)
endif(WIN32)

target_link_libraries(${CFD_TARGET_NAME} ${LINK_LIBRARY_LIST})

message(STATUS "Link library list: ${LINK_LIBRARY_LIST}.")

# For UGKS

#24. exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")

