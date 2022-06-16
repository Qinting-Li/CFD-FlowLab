cmake_minimum_required(VERSION 2.8)

#build the project:
project(PHengLEIv3d0 CXX)

#set working directory:
set(EXEC_WORKING_DIR PHengLEIv3d0)

#set COMPILER:
if(WIN32)
  if(CMAKE_GENERATOR MATCHES "Visual Studio 11 2012")
    SET(THIRD_PARTY_COMPILER_FLAG "libvs12")
    add_definitions(-DVS12)
  elseif(CMAKE_GENERATOR MATCHES "Visual Studio 12 2013")
    SET(THIRD_PARTY_COMPILER_FLAG "libvs12")
    add_definitions(-DVS12)
  elseif(CMAKE_GENERATOR MATCHES "Visual Studio 14 2015")
    SET(THIRD_PARTY_COMPILER_FLAG "libvs19")
  elseif(CMAKE_GENERATOR MATCHES "Visual Studio 15 2017")
    SET(THIRD_PARTY_COMPILER_FLAG "libvs19")
  elseif(CMAKE_GENERATOR MATCHES "Visual Studio 16 2019")
    SET(THIRD_PARTY_COMPILER_FLAG "libvs19")
  endif()
else()
  SET(CMAKE_C_COMPILER gcc)
  SET(THIRD_PARTY_COMPILER_FLAG "Unix")
endif()
include(${PROJECT_SOURCE_DIR}/../PHengLEI/PHMacro.cmake)
#set cmake environment:
set(CMAKE_USE_RELATIVE_PATHS ON)

#output message:
message(STATUS "Welcome to Build Project ${CMAKE_PROJECT_NAME}.")
message(STATUS "The OS is ${CMAKE_SYSTEM_ID}, version ${CMAKE_SYSTEM_VERSION}.")
message(STATUS "The COMPILER is ${THIRD_PARTY_COMPILER_FLAG}.")
message(STATUS "The working directory is ${CMAKE_BINARY_DIR}/${EXEC_WORKING_DIR}/${DLL_WORKING_DIR}.")
message(STATUS "--- Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

# Define some pre-definition.
opt(OCC "Use OpenCasCade Library" OFF)  # Use occ or not.
opt(DOUBLE_PRECISION "Use double precision or not" ON)
opt(SMARTARRAY "Use smart array more efficenet than origial array" ON)  # Use smart array or not.
opt_build(SpecDiffHybSolver "Build Spectra-Difference Hybrid Solver or not" OFF) #Build Spectra-Difference Hybrid Solver or not.
opt_build(DemoSolver "Build Solver Template or not" OFF)# Build Solver Template option.
opt_build(PBSolver "Build Pressure Based Solver or not" ON)# Build PBSolver option.
opt(OMP "Use OpenMP or not" OFF) # Use OPenMP or not.
opt(Windows_X64 "Compile PHengLEI X64 on Windows platform" ON)  ###
opt_build(API "Build PHengLEI framework or not" ON) #Build API or not.
opt_build(MESH "Build Mesh module or not" ON) #Build Mesh or not.
opt_build(CFD "Build CFD module or not" ON) #Build CFD or not.
opt(ErrorDataDebug "Build error data debug or not" OFF)

set(GPROF OFF)
 if (GPROF)
   set(CMAKE_CXX_FLAGS "$ENV{CMAKE_CXX_FLAGS} -pg")
   set(CMAKE_C_FLAGS "$ENV{CMAKE_C_FLAGS} -pg")
   set(CMAKE_EXE_LINKER_FLAGS "$ENV{CMAKE_EXE_LINKER_FLAGS} -pg")
   set(CMAKE_SHARED_LINKER_FLAGS "$ENV{CMAKE_SHARED_LINKER_FLAGS} -pg")
 endif()

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

if(BUILD_SpecDiffHybSolver)
  message(STATUS "BUILD_SpecDiffHybSolver is turned on!")
  add_definitions(-DUSE_SpecDiffHybSolver)
endif()

if(BUILD_DemoSolver)
  message(STATUS "BUILD_DemoSolver is turned on!")
  add_definitions(-DUSE_DEMOSOLVER)
endif()

if(BUILD_PBSolver)
  message(STATUS "BUILD_PBSolver is turned on!")
  add_definitions(-DUSE_PBSOLVER)
endif()

if(USE_OMP)
  message(STATUS "USE_OMP is turned on!")
  add_definitions(-DUSE_OMP)
  find_package (OpenMP)
  if (OpenMP_FOUND)
  message (STATUS "OpenMP_VERSION = ${OpenMP_VERSION}")
  message (STATUS "OpenMP_FOUND    = ${OpenMP_FOUND}")
  message (STATUS "OpenMP_CXX_FOUND  = ${OpenMP_CXX_FOUND}")
  message (STATUS "OpenMP_CXX_FLAGS  = ${OpenMP_CXX_FLAGS}")
  message (STATUS "OpenMP_CXX_INCLUDE_DIRS  = ${OpenMP_CXX_INCLUDE_DIRS}")
  message (STATUS "OpenMP_CXX_LIB_NAMES  = ${OpenMP_CXX_LIB_NAMES}")
  message (STATUS "OpenMP_CXX_LIBRARIES  = ${OpenMP_CXX_LIBRARIES}")
  message (STATUS "OpenMP_LIBRARY  = ${OpenMP_LIBRARY}")
  message (STATUS "OpenMP_CXX_SPEC_DATE  = ${OpenMP_CXX_SPEC_DATE}")
  message (STATUS "OpenMP_CXX_VERSION_MAJOR  = ${OpenMP_CXX_VERSION_MAJOR}")
  message (STATUS "OpenMP_CXX_VERSION_MINOR  = ${OpenMP_CXX_VERSION_MINOR}")
  message (STATUS "OpenMP_CXX_VERSION  = ${OpenMP_CXX_VERSION}") 
  endif (OpenMP_FOUND)
endif(USE_OMP)

if(USE_ErrorDataDebug)
  message(STATUS "USE_ErrorDataDebug is turned on!")
  add_definitions(-DUSE_ErrorDataDebug)
endif(USE_ErrorDataDebug)

if(WIN32)
if(USE_Windows_X64)
  message(STATUS "USE_Windows_X64 is turned on!")
  add_definitions(-DUSE_WINDOWS_X64)
endif(USE_Windows_X64)
endif(WIN32)

if(WIN32)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
endif()

# Add subdirectory or link library PHAPI.
set(API_TARGET_NAME "PHAPI")

includeHeadFiles()

set(HDF5_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/include)
include_directories("${HDF5_H_PATH}")
set(GMSH_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/Gmsh/include)
include_directories("${GMSH_H_PATH}")
set(PARMETIS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/parmetis/include)
include_directories("${PARMETIS_H_PATH}")
set(MEITIS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/metis/include)
include_directories("${MEITIS_H_PATH}")
set(CGNS_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/cgns/include)
include_directories("${CGNS_H_PATH}")

if(WIN32)
if(USE_Windows_X64)
  set(LIB_API_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  set(CMAKE_EXE_LINKER_FLAGS "/machine:x64")
else(USE_Windows_X64)
  set(LIB_API_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(USE_Windows_X64)
else(WIN32)
  set(LIB_API_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../API/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)

if(BUILD_API)
  message(STATUS "Build API!")
  add_subdirectory(../API ${EXEC_WORKING_DIR}/API)
if(WIN32)
  set_target_properties(${API_TARGET_NAME} PROPERTIES COMPILE_FLAGS "/MP")
endif(WIN32)
else(BUILD_API)
  message(STATUS "Link API!") 
  link_directories("${LIB_API_EXPORT_PATH}")  
endif(BUILD_API)

# Add subdirectory or link library PHMesh.
set(MESH_TARGET_NAME "PHMesh")
if(WIN32)
if(USE_Windows_X64)
  set(LIB_MESH_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../Mesh/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
else(USE_Windows_X64)
  set(LIB_MESH_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../Mesh/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(USE_Windows_X64)
else(WIN32)
  set(LIB_MESH_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../Mesh/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)

set(MESH_H_PATH ${PROJECT_SOURCE_DIR}/../Mesh/include)
include_directories(${MESH_H_PATH})
if(BUILD_MESH)
  message(STATUS "Build MESH!")
  add_subdirectory(../Mesh ${EXEC_WORKING_DIR}/Mesh)
if(WIN32)
  set_target_properties(${Mesh_TARGET_NAME} PROPERTIES COMPILE_FLAGS "/MP")
endif(WIN32)
else(BUILD_MESH)
  message(STATUS "Link MESH!") 
  link_directories("${LIB_MESH_EXPORT_PATH}")  
endif(BUILD_MESH)


# Add subdirectory or link library PHCFD.
set(CFD_TARGET_NAME "PHCFD")
if(WIN32)
if(USE_Windows_X64)
  set(LIB_CFD_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
else(USE_Windows_X64)
  set(LIB_CFD_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(USE_Windows_X64)
else(WIN32)
  set(LIB_CFD_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)

set(CFD_H_PATH   ${PROJECT_SOURCE_DIR}/../CFD)
include_directories(${CFD_H_PATH}/NSSolver/include)
include_directories(${CFD_H_PATH}/Solver/include)
include_directories(${CFD_H_PATH}/TurbSolver/include)
include_directories(${CFD_H_PATH}/PostProcess/include)
include_directories(${CFD_H_PATH}/HOSolverUnstruct/include)

if(BUILD_SpecDiffHybSolver)
include_directories(${CFD_H_PATH}/SpecDiffHybSolver/include)
endif(BUILD_SpecDiffHybSolver)

if(BUILD_DemoSolver)
include_directories(${CFD_H_PATH}/DemoSolver/include)
endif(BUILD_DemoSolver)

if(BUILD_PBSolver)
include_directories(${CFD_H_PATH}/PBSolver/include)
endif(BUILD_PBSolver)
if(BUILD_CFD)
  message(STATUS "Build CFD!")
  add_subdirectory(../CFD ${EXEC_WORKING_DIR}/CFD)
if(WIN32)
  set_target_properties(${CFD_TARGET_NAME} PROPERTIES COMPILE_FLAGS "/MP")
endif(WIN32)
else(BUILD_CFD)
  message(STATUS "Link CFD!") 
  link_directories("${LIB_CFD_EXPORT_PATH}")  
endif(BUILD_CFD)

if(WIN32)
if(USE_Windows_X64)
  link_directories("${PROJECT_SOURCE_DIR}/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x64")
else(USE_Windows_X64)
  link_directories("${PROJECT_SOURCE_DIR}/3rdparty/${THIRD_PARTY_COMPILER_FLAG}/win-x86")
endif(USE_Windows_X64)
message(STATUS "CMAKE_CXX_COMPILER! ${CMAKE_CXX_COMPILER}")
else(WIN32)
  link_directories("${PROJECT_SOURCE_DIR}/3rdparty/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
endif(WIN32)

#look for CMakeLists.txt in subdirectory <codes> and build an executable target:
add_subdirectory(../Main ${EXEC_WORKING_DIR})

#exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
message(STATUS "Finished build project ${PROJECT_NAME}!")


