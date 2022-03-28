cmake_minimum_required(VERSION 2.8)

#build the project:
project (KDTree CXX) 
 
include(${PROJECT_SOURCE_DIR}/../../PHMacro.cmake)

add_definitions(-DPH_PARALLEL )
#set working directory:
#set(EXEC_WORKING_DIR KDTree)

#set cmake environment:
#set(CMAKE_USE_RELATIVE_PATHS ON)
opt(Windows_X64 "Compile PHengLEI X64 on Windows platform" OFF)  ###
#output message:
#message(STATUS "Welcome to Build Project ${CMAKE_PROJECT_NAME}.")
#message(STATUS "The OS is ${CMAKE_SYSTEM_NAME}, version ${CMAKE_SYSTEM_VERSION}.")
#message(STATUS "The working directory is ${CMAKE_BINARY_DIR}/${EXEC_WORKING_DIR}/${DLL_WORKING_DIR}.")
#message(STATUS "--- Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
#message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
#message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

## Add subdirectory or link library.
set(DataStruct_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/DataStruct/include)
set(PHMPI_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/FYMPI/include)
include_directories("${DataStruct_H_PATH}")
include_directories("${PHMPI_H_PATH}")

#3. set environment for MPICH:
set(THIRDPARTY_LIB_DIR)
if(WIN32)
  #set(MPICH_HOME_DIRECTORY "C:/Program Files (x86)/Microsoft SDKs/MPI/Include" CACHE PATH "path to installed MPICH")
  #link_directories("${MPICH_HOME_DIRECTORY}/lib")
  #include_directories("${MPICH_HOME_DIRECTORY}/include")
  set(MSMPI_HOME_DIRECTORY "C:/Program Files (x86)/Microsoft SDKs/MPI" CACHE PATH "path to installed MPICH")
  include_directories("${MSMPI_HOME_DIRECTORY}/Include") 
if(USE_Windows_X64)
  set(LIB_DataStruct_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/DataStruct/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_FYMPI_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/FYMPI/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 
else(USE_Windows_X64)
  set(LIB_DataStruct_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/DataStruct/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
   set(LIB_FYMPI_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/FYMPI/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(USE_Windows_X64)
else(WIN32)
  set(CMAKE_CXX_COMPILER  mpicxx)
  set(CMAKE_C_COMPILER    mpicc)
  set(LIB_DataStruct_EXPORT_PATH  ${PROJECT_SOURCE_DIR}/../../../API0/DataStruct/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_FYMPI_EXPORT_PATH       ${PROJECT_SOURCE_DIR}/../../../API0/FYMPI/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)


message(STATUS "The DataStruct path is: ${LIB_DataStruct_EXPORT_PATH}.") 
link_directories("${LIB_DataStruct_EXPORT_PATH}")

message(STATUS "The FYMPI path is: ${LIB_FYMPI_EXPORT_PATH}.") 
link_directories("${LIB_FYMPI_EXPORT_PATH}")
add_definitions(-DFY_PARALLEL -DFY_CGNS -DFY_DOUBLE_PREC)

set(DataStruct_TARGET_NAME "PHDataStruct")

include_directories(${DataStruct_H_PATH})

list(APPEND LINK_LIBRARY_LIST  PHDataStruct PHMPI)

add_executable(KDTree src/KDTreeTest.cpp)

#23. link a target to given libraries:
target_link_libraries(KDTree ${LINK_LIBRARY_LIST})

message(STATUS " --- Done ! ---")
