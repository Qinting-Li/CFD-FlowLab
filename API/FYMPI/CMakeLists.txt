#! @file      CMakeLists.txt
#! @brief     It defines the PHMPI.
#! @author    He 
cmake_minimum_required(VERSION 2.6)

#1. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory <FYMPI>, Start to call ${CMAKE_CURRENT_LIST_FILE} ---")
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

set(PHMPI_TARGET_NAME "PHMPI")
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
	set (BUILD_LIBRARY_TYPE "Static" CACHE STRING "The type of PHMPI libraries" FORCE)
  SET_PROPERTY(CACHE BUILD_LIBRARY_TYPE PROPERTY STRINGS Shared Static)
endif()

#3. set environment for MPICH and 3rdParty lib directory:
if(WIN32)
  set(MSMPI_HOME_DIRECTORY "C:/Program Files (x86)/Microsoft SDKs/MPI" CACHE PATH "path to installed MPICH")
  include_directories("${MSMPI_HOME_DIRECTORY}/Include")
  if(USE_Windows_X64)
	message(STATUS "USE_Windows_X64 is turned on!")
	add_definitions(-DUSE_WINDOWS_X64)
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
${MAIN_H_PATH}/PHMpi.h
${MAIN_H_PATH}/Parallel.h
${MAIN_H_PATH}/ActionKey.h
)
source_group(include FILES ${MAIN_H_LIST}) 
include_directories(${MAIN_H_PATH})

#7. set directory <src>:
set(MAIN_CPP_PATH  src)
set(MAIN_CPP_LIST
${MAIN_CPP_PATH}/PHMpi.cpp
${MAIN_CPP_PATH}/Parallel.cpp
${MAIN_CPP_PATH}/ActionKey.cpp
)
source_group(src FILES ${MAIN_CPP_LIST}) 
include_directories(${MAIN_CPP_PATH})

#20. set header file list:
set(HEAD_LIST ${MAIN_H_LIST})

#21. set source file list:
set(CPP_LIST ${MAIN_CPP_LIST})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../API/Common/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../API/DataStruct/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../API/Data/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../API/Math/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../API/IO/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../../API/Toolkit/include")
add_definitions(-DPH_PARALLEL -DPH_CGNS)

#22. build executable program or library to the project using the specified source files and header files:
if ("${BUILD_LIBRARY_TYPE}" STREQUAL "Shared")
  # Building Shared Library.   
  message(STATUS "Building Shared Library ${PHMPI_TARGET_NAME} ...")
  add_library(${PHMPI_TARGET_NAME} SHARED ${CPP_LIST} ${HEAD_LIST})
elseif ("${BUILD_LIBRARY_TYPE}" STREQUAL "Static")
  # Building Static Library.
  message(STATUS "Building Static Library ${PHMPI_TARGET_NAME} ...")
  add_library(${PHMPI_TARGET_NAME} STATIC ${CPP_LIST} ${HEAD_LIST}) 
endif() 
# Collect the linraries to be linked.
list(APPEND LINK_LIBRARY_LIST DataStruct)
if(WIN32)
  list(APPEND LINK_LIBRARY_LIST msmpi)
endif(WIN32)

message(STATUS "Link library list: ${LINK_LIBRARY_LIST}.")

#23. link a target to given libraries:
target_link_libraries(${PHMPI_TARGET_NAME} )
	
#24. exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
