#! @file      CMakeLists.txt
#! @brief     It defines the Math.
#! @author    He Xianyao.
cmake_minimum_required(VERSION 2.6)

#1. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory <Math>, Start to call ${CMAKE_CURRENT_LIST_FILE} ---")
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

# Define some pre-definition.
include(PHMacro.cmake)
if(USE_SMARTARRAY)
  message(STATUS "USE_SMARTARRAY is turned on!")
  add_definitions(-DUSE_SMARTARRAY)
endif()
if(USE_DOUBLE_PRECISION)
  message(STATUS "USE_DOUBLE_PRECISION is turned on!")
  add_definitions(-DPH_DOUBLE_PREC)
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

set(Math_TARGET_NAME "PHMath")
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
if (NOT BUILD_LIBRARY_TYPE)
	set (BUILD_LIBRARY_TYPE "Static" CACHE STRING "The type of Math libraries" FORCE)
  SET_PROPERTY(CACHE BUILD_LIBRARY_TYPE PROPERTY STRINGS Shared Static)
endif()

#3. set environment for MPICH:
if(NOT WIN32)
  set(CMAKE_CXX_COMPILER  mpicxx)
  set(CMAKE_C_COMPILER    mpicc)
endif()

#6. set subdirectory <include>:
set(MAIN_H_PATH  include)
set(MAIN_H_LIST
${MAIN_H_PATH}/Math_BasisFunction.h
${MAIN_H_PATH}/Math_Limiter.h
)
source_group(include FILES ${MAIN_H_LIST}) 
include_directories(${MAIN_H_PATH})

#7. set directory <src>:
set(MAIN_CPP_PATH  src)
set(MAIN_CPP_LIST
${MAIN_CPP_PATH}/Math_BasisFunction.cpp
)
source_group(src FILES ${MAIN_CPP_LIST}) 
include_directories(${MAIN_CPP_PATH})

#20. set header file list:
set(HEAD_LIST ${MAIN_H_LIST})

#21. set source file list:
set(CPP_LIST ${MAIN_CPP_LIST})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../API/Common/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../API/DataStruct/include")
add_definitions(-DPH_PARALLEL -DPH_CGNS)

#22. build executable program or library to the project using the specified source files and header files:
if ("${BUILD_LIBRARY_TYPE}" STREQUAL "Shared")
  # Building Shared Library.   
  message(STATUS "Building Shared Library ${Math_TARGET_NAME} ...")
  add_library(${Math_TARGET_NAME} SHARED ${CPP_LIST} ${HEAD_LIST})
elseif ("${BUILD_LIBRARY_TYPE}" STREQUAL "Static")
  # Building Static Library.
  message(STATUS "Building Static Library ${Math_TARGET_NAME} ...")
  add_library(${Math_TARGET_NAME} STATIC ${CPP_LIST} ${HEAD_LIST}) 
endif() 

#23. link a target to given libraries:
target_link_libraries(${Math_TARGET_NAME} )
	
#24. exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
