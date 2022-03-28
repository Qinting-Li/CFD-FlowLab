cmake_minimum_required(VERSION 2.8)

#build the project:
project (GeoNodeTopo_str CXX) 
 
include(${PROJECT_SOURCE_DIR}/../../PHMacro.cmake)

add_definitions(-DPH_PARALLEL )

#set working directory:
set(EXEC_WORKING_DIR GeoNodeTopo_str)

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
set(Data_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Data/include)
set(Math_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Math/include)
set(PHMPI_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/FYMPI/include)
set(Toolkit_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Toolkit/include)
set(IO_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/IO/include)
set(Geometry_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Geometry/include)
set(PreProcess_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/PreProcess/include)
set(PostProcess_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/PostProcess/include)
set(INCLUDE_CFD_PATH ${PROJECT_SOURCE_DIR}/../../../CFD/main)
set(MESH_H_PATH   ${PROJECT_SOURCE_DIR}/../../../Mesh/include)
set(Common_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Common/include)
set(Physics_H_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Physics/include)
set(CGNS_H_PATH ${PROJECT_SOURCE_DIR}/../../../PHengLEI/3rdparty/cgns/include)
set(HDF5_H_PATH ${PROJECT_SOURCE_DIR}/../../../PHengLEI/3rdparty/HDF5/include)
include_directories("${DataStruct_H_PATH}")
include_directories("${Data_H_PATH}")
include_directories("${Math_H_PATH}")
include_directories("${PHMPI_H_PATH}")
include_directories("${Toolkit_H_PATH}")
include_directories("${IO_H_PATH}")
include_directories("${Geometry_H_PATH}")
include_directories("${PreProcess_H_PATH}")
include_directories("${PostProcess_H_PATH}")
include_directories("${MESH_H_PATH}")
include_directories("${Common_H_PATH}")
include_directories("${Physics_H_PATH}")
include_directories("${CGNS_H_PATH}")
include_directories("${HDF5_H_PATH}")
message(STATUS "----HDF5 ${HDF5_H_PATH}   --- ")

#3. set environment for MPICH:
set(THIRDPARTY_LIB_DIR)
if(WIN32)
  #set(MPICH_HOME_DIRECTORY "C:/Program Files (x86)/Microsoft SDKs/MPI/Include" CACHE PATH "path to installed MPICH")
  #link_directories("${MPICH_HOME_DIRECTORY}/lib")
  #include_directories("${MPICH_HOME_DIRECTORY}/include")
  set(MSMPI_HOME_DIRECTORY "C:/Program Files (x86)/Microsoft SDKs/MPI" CACHE PATH "path to installed MPICH")
  include_directories("${MSMPI_HOME_DIRECTORY}/Include") 
if(USE_Windows_X64)
  set(LIB_API_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_DataStruct_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/DataStruct/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  #  set(LIB_Data_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Data/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_Math_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Math/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_FYMPI_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/FYMPI/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_Toolkit_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Toolkit/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_IO_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/IO/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(CMAKE_EXE_LINKER_FLAGS "/machine:x64")
  link_directories("${MSMPI_HOME_DIRECTORY}/Lib/x64")
  set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/../../3rdparty/win-x64")
  set(LIB_MESH_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../Mesh/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})


  set(LIB_DATA_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Data/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_GEOMETRY_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Geometry/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_COMMON_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Common/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_POSTPROCESS_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/PostProcess/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_PREPROCESS_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/PreProcess/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_PHYSICS_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Physics/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(LIB_CFD_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../../API0/CFD/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 


else(USE_Windows_X64)
   set(LIB_API_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
   set(LIB_DataStruct_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/DataStruct/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
   #   set(LIB_Data_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Data/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
   set(LIB_Math_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Math/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
   set(LIB_FYMPI_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/FYMPI/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
   set(LIB_Toolkit_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Toolkit/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
   set(LIB_IO_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/IO/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
   link_directories("${MSMPI_HOME_DIRECTORY}/Lib/x86")
   set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/../../3rdparty/win-x86")
   set(LIB_MESH_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../Mesh/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})  

  set(LIB_DATA_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Data/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_GEOMETRY_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Geometry/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_COMMON_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Common/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_POSTPROCESS_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/PostProcess/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_PREPROCESS_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/PreProcess/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_PHYSICS_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/Physics/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  set(LIB_CFD_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../../API0/CFD/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})   
endif(USE_Windows_X64)
else(WIN32)
  set(CMAKE_CXX_COMPILER  mpicxx)
  set(CMAKE_C_COMPILER    mpicc)
  set(LIB_DataStruct_EXPORT_PATH  ${PROJECT_SOURCE_DIR}/../../../API0/DataStruct/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_Math_EXPORT_PATH        ${PROJECT_SOURCE_DIR}/../../../API0/Math/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_FYMPI_EXPORT_PATH       ${PROJECT_SOURCE_DIR}/../../../API0/FYMPI/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_Toolkit_EXPORT_PATH     ${PROJECT_SOURCE_DIR}/../../../API0/Toolkit/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(THIRDPARTY_LIB_DIR "${PROJECT_SOURCE_DIR}/../../../PHengLEI/3rdparty/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
  set(LIB_MESH_EXPORT_PATH        ${PROJECT_SOURCE_DIR}/../../../Mesh/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  
  set(LIB_DATA_EXPORT_PATH        ${PROJECT_SOURCE_DIR}/../../../API0/Data/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_IO_EXPORT_PATH          ${PROJECT_SOURCE_DIR}/../../../API0/IO/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_GEOMETRY_EXPORT_PATH    ${PROJECT_SOURCE_DIR}/../../../API0/Geometry/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_COMMON_EXPORT_PATH      ${PROJECT_SOURCE_DIR}/../../../API0/Common/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_POSTPROCESS_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../../../API0/PostProcess/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_PREPROCESS_EXPORT_PATH  ${PROJECT_SOURCE_DIR}/../../../API0/PreProcess/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_PHYSICS_EXPORT_PATH     ${PROJECT_SOURCE_DIR}/../../../API0/Physics/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  set(LIB_CFD_EXPORT_PATH         ${PROJECT_SOURCE_DIR}/../../../CFD/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)

message(STATUS "The DataStruct path is: ${LIB_DataStruct_EXPORT_PATH}.") 
link_directories("${LIB_DataStruct_EXPORT_PATH}")

message(STATUS "The Math path is: ${LIB_Math_EXPORT_PATH}.") 
link_directories("${LIB_Math_EXPORT_PATH}")

message(STATUS "The FYMPI path is: ${LIB_FYMPI_EXPORT_PATH}.") 
link_directories("${LIB_FYMPI_EXPORT_PATH}")

message(STATUS "The Toolkit path is: ${LIB_Toolkit_EXPORT_PATH}.") 
link_directories("${LIB_Toolkit_EXPORT_PATH}")

message(STATUS "The 3rdparty path is: ${THIRDPARTY_LIB_DIR}.") 
link_directories("${THIRDPARTY_LIB_DIR}")

message(STATUS "The mesh path is: ${LIB_MESH_EXPORT_PATH}.") 
link_directories("${LIB_MESH_EXPORT_PATH}")  

message(STATUS "The DATA path is: ${LIB_DATA_EXPORT_PATH}.") 
link_directories("${LIB_DATA_EXPORT_PATH}")

message(STATUS "The IO path is: ${LIB_IO_EXPORT_PATH}.") 
link_directories("${LIB_IO_EXPORT_PATH}")

message(STATUS "The Geometry path is: ${LIB_GEOMETRY_EXPORT_PATH}.") 
link_directories("${LIB_GEOMETRY_EXPORT_PATH}")

message(STATUS "The Common path is: ${PROJECT_SOURCE_DIR}/../../../API0/Common/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}.") 
link_directories("${PROJECT_SOURCE_DIR}/../../../API0/Common/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")

message(STATUS "The PostProcess path is: ${PROJECT_SOURCE_DIR}/../../../API0/PostProcess/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}.") 
link_directories("${PROJECT_SOURCE_DIR}/../../../API0/PostProcess/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")

message(STATUS "The PreProcess path is: ${PROJECT_SOURCE_DIR}/../../../API0/PreProcess/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}.") 
link_directories("${PROJECT_SOURCE_DIR}/../../../API0/PreProcess/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")

message(STATUS "The PHYSICS path is: ${LIB_PHYSICS_EXPORT_PATH}.") 
link_directories("${LIB_PHYSICS_EXPORT_PATH}") 

message(STATUS "The CFD path is: ${LIB_CFD_EXPORT_PATH}.") 
link_directories("${LIB_CFD_EXPORT_PATH}") 

add_definitions(-DFY_PARALLEL -DFY_CGNS -DFY_DOUBLE_PREC)

#link library
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

message(STATUS "Linking ${MESH_TARGET_NAME} ...")
list(APPEND LINK_LIBRARY_LIST ${MESH_TARGET_NAME})

message(STATUS "Linking ${CFD_TARGET_NAME} ...")
list(APPEND LINK_LIBRARY_LIST ${CFD_TARGET_NAME})

list(APPEND LINK_LIBRARY_LIST ${Physics_TARGET_NAME} ${PreProcess_TARGET_NAME} ${PostProcess_TARGET_NAME} ${Common_TARGET_NAME} ${Geometry_TARGET_NAME} ${IO_TARGET_NAME} ${Toolkit_TARGET_NAME} ${PHMPI_TARGET_NAME} ${DataStruct_TARGET_NAME} ${Data_TARGET_NAME} ${Math_TARGET_NAME} PHPreProcess PHPhysics PHMPI PHToolkit metis parmetis cgns hdf5 tecio)
if(WIN32)
  list(APPEND LINK_LIBRARY_LIST msmpi)
endif(WIN32)

add_executable(UT_Actionkey src/UT_Actionkey.cpp)

#23. link a target to given libraries:
target_link_libraries(UT_Actionkey ${LINK_LIBRARY_LIST})

message(STATUS " --- Done ! ---")
