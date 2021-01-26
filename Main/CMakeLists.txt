#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +  
#          P   P  H   H  E      NN   N  G      L      E       I          +
#          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
#          P      H   H  E      N  N N  G   G  L      E       I          +
#          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
#------------------------------------------------------------------------+
#          Platform for Hybrid Engineering Simulation of Flows           +
#          China Aerodynamics Research and Development Center            +
#                     (C) Copyright, Since 2010                          +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! @file      CMakeLists.txt
#! @brief     It defines the project.
#! @author    Bell, Xu Qingxin, He Kun, He Xin.
cmake_minimum_required(VERSION 2.8)

#1. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory <Main>, Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

#2. set executable name:
set(EXECUTABLE_TARGET_NAME "${CMAKE_PROJECT_NAME}")
message(STATUS "The executable name is ${EXECUTABLE_TARGET_NAME}.")

# Define some pre-definition.
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

#5. set current directory:

#6. set subdirectory <include>:
set(MAIN_H_PATH   include)
set(MAIN_H_LIST
${MAIN_H_PATH}/CFDSimulation.h
${MAIN_H_PATH}/Simulation.h
${MAIN_H_PATH}/GridFactory.h
${MAIN_H_PATH}/FieldViewIn.h
${MAIN_H_PATH}/GridPartition.h
${MAIN_H_PATH}/GridView.h
${MAIN_H_PATH}/FluentIn.h
${MAIN_H_PATH}/Gmesh.h
${MAIN_H_PATH}/GridgenIn.h
${MAIN_H_PATH}/SolversList.h
)
source_group(Main\\include FILES ${MAIN_H_LIST}) 
include_directories(${MAIN_H_PATH})

#7. set directory <src>:
set(MAIN_CPP_PATH  src)
set(MAIN_CPP_LIST
${MAIN_CPP_PATH}/main.cpp
${MAIN_CPP_PATH}/CFDSimulation.cpp
${MAIN_CPP_PATH}/Simulation.cpp
${MAIN_CPP_PATH}/GridFactory.cpp
${MAIN_CPP_PATH}/FieldViewIn.cpp
${MAIN_CPP_PATH}/GridPartition.cpp
${MAIN_CPP_PATH}/GridView.cpp
${MAIN_CPP_PATH}/FluentIn.cpp
${MAIN_CPP_PATH}/Gmesh.cpp
${MAIN_CPP_PATH}/GridgenIn.cpp
)

source_group(Main\\src FILES ${MAIN_CPP_LIST}) 
include_directories(${MAIN_CPP_PATH})

#20.1 set header file list
set(HEAD_LIST 
${MAIN_H_LIST}
)

#20.2 set source file list:
set(CPP_LIST 
${MAIN_CPP_LIST}
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
endif(USE_OCC)

#22. build executable program:

# Building Executable Program.
message(STATUS "Building Executable Program ${EXECUTABLE_TARGET_NAME} ...")
add_executable(${EXECUTABLE_TARGET_NAME} ${CPP_LIST} ${HEAD_LIST})  
if(WIN32) 
    set_target_properties(${EXECUTABLE_TARGET_NAME} PROPERTIES COMPILE_FLAGS "/MP")
endif(WIN32)
# Collect the linraries to be linked.
# link the API/MESH etc. libraries, no matter build them or not.
# Important: the top library should be list firstly!!!
# For example, the Mesh is previous API.

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
set(MESH_TARGET_NAME "PHMesh")
set(CFD_TARGET_NAME "PHCFD")

if(BUILD_LBM)
  set(LBM_TARGET_NAME "PHLBM")
endif(BUILD_LBM)


message(STATUS "Linking ${CFD_TARGET_NAME} ...")
list(APPEND LINK_LIBRARY_LIST ${CFD_TARGET_NAME})

list(APPEND LINK_LIBRARY_LIST ${Physics_TARGET_NAME} ${PreProcess_TARGET_NAME} ${PostProcess_TARGET_NAME} ${Common_TARGET_NAME} ${Geometry_TARGET_NAME} ${IO_TARGET_NAME} ${Toolkit_TARGET_NAME} ${PHMPI_TARGET_NAME} ${DataStruct_TARGET_NAME} ${Data_TARGET_NAME} ${Math_TARGET_NAME})

message(STATUS "Link library list: ${LINK_LIBRARY_LIST}.")


#23. link a target to given libraries:
target_link_libraries(${EXECUTABLE_TARGET_NAME} ${LINK_LIBRARY_LIST})


# Some pre-definitions that occ need.
if(USE_OCC)
  target_link_libraries(${EXECUTABLE_TARGET_NAME} ${OCC_LIBS})  
  add_definitions(-DWNT)
  add_definitions(-DNOTCL)
endif(USE_OCC)
	
# For UGKS

#24. exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")

