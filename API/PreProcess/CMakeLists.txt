#! @file      CMakeLists.txt
#! @brief     It defines the PreProcess.
#! @author    He Xianyao.
cmake_minimum_required(VERSION 2.6)

#1. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory <PreProcess>, Start to call ${CMAKE_CURRENT_LIST_FILE} ---")
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

set(PreProcess_TARGET_NAME "PHPreProcess")
if(WIN32)
if(USE_Windows_X64)
	set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(CMAKE_EXE_LINKER_FLAGS "/machine:x64")
else(USE_Windows_X64)
	set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(USE_Windows_X64)
else(WIN32)
	set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)
#set options
if(NOT BUILD_LIBRARY_TYPE)
	set (BUILD_LIBRARY_TYPE "Static" CACHE STRING "The type of PreProcess libraries" FORCE)
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

#6. set subdirectory <include>:
set(MAIN_H_PATH  include)
set(MAIN_H_LIST
${MAIN_H_PATH}/Pre_Block_Struct.h
${MAIN_H_PATH}/Pre_UnstructGridPartition.h
${MAIN_H_PATH}/Pre_WalldistCompute.h
${MAIN_H_PATH}/Pre_GridConversion.h
${MAIN_H_PATH}/Pre_FieldViewToPHengLEI_Unstruct.h
${MAIN_H_PATH}/Pre_CGNSToPHengLEI_Unstruct.h
${MAIN_H_PATH}/Pre_CGNSToPlot3D_Struct.h
${MAIN_H_PATH}/Pre_FluentToPHengLEI_Unstruct.h
${MAIN_H_PATH}/Pre_HDF5File.h
${MAIN_H_PATH}/WallDist.h
${MAIN_H_PATH}/Pre_GridMirror.h
${MAIN_H_PATH}/Pre_GridRefine.h
${MAIN_H_PATH}/Pre_StrGridPartition.h
${MAIN_H_PATH}/GridTranslator.h
${MAIN_H_PATH}/CombinGrid.h
${MAIN_H_PATH}/GridBase.h
${MAIN_H_PATH}/MixGrid.h
)
source_group(include FILES ${MAIN_H_LIST}) 
include_directories(${MAIN_H_PATH})

#7. set directory <src>:
set(MAIN_CPP_PATH  src)
set(MAIN_CPP_LIST
${MAIN_CPP_PATH}/Pre_Block_Struct.cpp
${MAIN_CPP_PATH}/Pre_UnstructGridPartition.cpp
${MAIN_CPP_PATH}/Pre_WalldistCompute.cpp
${MAIN_CPP_PATH}/Pre_GridConversion.cpp
${MAIN_CPP_PATH}/Pre_FieldViewToPHengLEI_Unstruct.cpp
${MAIN_CPP_PATH}/Pre_CGNSToPHengLEI_Unstruct.cpp
${MAIN_CPP_PATH}/Pre_CGNSToPlot3D_Struct.cpp
${MAIN_CPP_PATH}/Pre_FluentToPHengLEI_Unstruct.cpp
${MAIN_CPP_PATH}/Pre_HDF5File.cpp
${MAIN_CPP_PATH}/WallDist.cpp
${MAIN_CPP_PATH}/Pre_GridMirror.cpp
${MAIN_CPP_PATH}/Pre_GridRefine.cpp
${MAIN_CPP_PATH}/Pre_StrGridPartition.cpp
${MAIN_CPP_PATH}/GridTranslator.cpp
${MAIN_CPP_PATH}/CombinGrid.cpp
${MAIN_CPP_PATH}/GridBase.cpp
${MAIN_CPP_PATH}/MixGrid.cpp
)
source_group(src FILES ${MAIN_CPP_LIST}) 
include_directories(${MAIN_CPP_PATH})

#20. set header file list:
set(HEAD_LIST ${MAIN_H_LIST})

#20. set source file list:
set(CPP_LIST ${MAIN_CPP_LIST})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../Common/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../Geometry/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../Physics/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../Data/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../DataStruct/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../Math/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../FYMPI/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../Toolkit/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../IO/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../PostProcess/include")
set(HDF5_H_PATH      ${CMAKE_CURRENT_SOURCE_DIR}/../../PHengLEI/3rdparty/HDF5/include)
include_directories("${HDF5_H_PATH}")

if(USE_TecplotLib)
set(TECIO_H_PATH     ${CMAKE_CURRENT_SOURCE_DIR}/../../PHengLEI/3rdparty/TECIO/include)
include_directories("${TECIO_H_PATH}")
endif(USE_TecplotLib)

set(PARMETIS_H_PATH  ${CMAKE_CURRENT_SOURCE_DIR}/../../PHengLEI/3rdparty/parmetis/include)
include_directories("${PARMETIS_H_PATH}")
set(MEITIS_H_PATH    ${CMAKE_CURRENT_SOURCE_DIR}/../../PHengLEI/3rdparty/metis/include)
include_directories("${MEITIS_H_PATH}")
set(CGNS_H_PATH      ${CMAKE_CURRENT_SOURCE_DIR}/../../PHengLEI/3rdparty/cgns/include)
include_directories("${CGNS_H_PATH}")
add_definitions(-DPH_PARALLEL -DPH_CGNS)

#22. build executable program or library to the project using the specified source files and header files:
if ("${BUILD_LIBRARY_TYPE}" STREQUAL "Shared")
  # Building Shared Library.   
  message(STATUS "Building Shared Library ${PreProcess_TARGET_NAME} ...")
  add_library(${PreProcess_TARGET_NAME} SHARED ${CPP_LIST} ${HEAD_LIST})
elseif ("${BUILD_LIBRARY_TYPE}" STREQUAL "Static")
  # Building Static Library.
  message(STATUS "Building Static Library ${PreProcess_TARGET_NAME} ...")
  add_library(${PreProcess_TARGET_NAME} STATIC ${CPP_LIST} ${HEAD_LIST}) 
endif() 

#set(Physics_TARGET_NAME "PHPhysics")
list(APPEND LINK_LIBRARY_LIST metis parmetis cgns hdf5)

if(USE_TecplotLib)
  list(APPEND LINK_LIBRARY_LIST tecio)
endif(USE_TecplotLib)

if(WIN32)
  list(APPEND LINK_LIBRARY_LIST msmpi)
endif(WIN32)

#23. link a target to given libraries:
target_link_libraries(${PreProcess_TARGET_NAME} )
	
#24. exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
