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
  add_compile_options(-std=c++14)
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
opt(ARRAY_OPT "Use array optimization for structure solver" OFF)  # Use array optimization or not.
opt(SpecDiffHybSolver "Use Spectra-Difference Hybrid Solver or not" OFF) #Use Spectra-Difference Hybrid Solver or not.
opt(NSDEMO "Use Solver Template or not" OFF) # Use Solver Template option
opt(OMP "Use OpenMP or not" OFF) # Use OPenMP or not.
opt(Windows_X64 "Compile PHengLEI X64 on Windows platform" ON)  ###
opt_build(API "Build PHengLEI framework or not" ON) #Build API or not.
opt_build(MESH "Build Mesh module or not" OFF) #Build Mesh or not.
opt_build(CFD "Build CFD module or not" ON) #Build CFD or not.
opt_build(TUTORIAL "Build Tutorial or not" OFF) #Build API tutorials or not.
opt(TecplotLib "Use tecplot lib to dump binary flow file" OFF) #Use tecplotlib or not.
opt_build(LBM "Build LBM module or not" OFF) #Build LBM or not.
opt(ErrorDataDebug "Build error data debug or not" OFF)

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
if(USE_SpecDiffHybSolver)
  message(STATUS "USE_SpecDiffHybSolver is turned on!")
  add_definitions(-DUSE_SpecDiffHybSolver)
endif()

if(USE_NSDEMO)
  message(STATUS "USE_NSDEMO is turned on!")
  add_definitions(-DUSE_NSDEMO)
endif()

if(USE_OMP)
  message(STATUS "USE_OMP is turned on!")
  add_definitions(-DUSE_OMP)
endif()

if(USE_ErrorDataDebug)
  message(STATUS "USE_ErrorDataDebug is turned on!")
  add_definitions(-DUSE_ErrorDataDebug)
endif(USE_ErrorDataDebug)


if(USE_TecplotLib)
  message(STATUS "USE_TecplotLib is turned on!")
  add_definitions(-DUSE_TecplotLib)
endif(USE_TecplotLib)

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
if(USE_TecplotLib)
  set(TECIO_H_PATH ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/TECIO/include)
  include_directories("${TECIO_H_PATH}")
endif(USE_TecplotLib)
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


if(BUILD_LBM) 
  # Add subdirectory or link library PHLBM.
  set(LBM_TARGET_NAME "PHLBM")
  if(WIN32)
  if(USE_Windows_X64)
    set(LIB_LBM_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../LBM/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
  else(USE_Windows_X64)
    set(LIB_LBM_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../LBM/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
  endif(USE_Windows_X64)
  else(WIN32)
    set(LIB_LBM_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../LBM/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
  endif(WIN32)
  
  set(LBM_H_PATH   ${PROJECT_SOURCE_DIR}/../LBM/include)
  include_directories(${LBM_H_PATH})
  if(BUILD_LBM)
    message(STATUS "Build LBM!")
    add_subdirectory(../LBM ${EXEC_WORKING_DIR}/LBM)
    if(WIN32)
      set_target_properties(${LBM_TARGET_NAME} PROPERTIES COMPILE_FLAGS "/MP")
    endif(WIN32)
  else(BUILD_LBM)
    message(STATUS "Link LBM!") 
    link_directories("${LIB_LBM_EXPORT_PATH}")  
  endif(BUILD_LBM)
endif(BUILD_LBM) 

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
if(USE_SpecDiffHybSolver)
include_directories(${CFD_H_PATH}/SpecDiffHybSolver/include)
endif(USE_SpecDiffHybSolver)
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

if(BUILD_TUTORIAL)
  ADD_SUBDIRECTORY(Tutorial/CGNSConversion            ${EXEC_WORKING_DIR}/Tutorial/CGNSConversion           )
  ADD_SUBDIRECTORY(Tutorial/GridReader                ${EXEC_WORKING_DIR}/Tutorial/GridReader               )
  ADD_SUBDIRECTORY(Tutorial/MathBasicTest             ${EXEC_WORKING_DIR}/Tutorial/MathBasicTest            )
  ADD_SUBDIRECTORY(Tutorial/LimiterFuncTest           ${EXEC_WORKING_DIR}/Tutorial/LimiterFuncTest          )
  ADD_SUBDIRECTORY(Tutorial/BasisFuncTest             ${EXEC_WORKING_DIR}/Tutorial/BasisFuncTest            )
  ADD_SUBDIRECTORY(Tutorial/FieldViewGridConversion   ${EXEC_WORKING_DIR}/Tutorial/FieldViewGridConversion  )
  ADD_SUBDIRECTORY(Tutorial/PreUnstructGridPartition  ${EXEC_WORKING_DIR}/Tutorial/PreUnstructGridPartition )
  ADD_SUBDIRECTORY(Tutorial/PreStructGridPartition    ${EXEC_WORKING_DIR}/Tutorial/PreStructGridPartition   )
  ADD_SUBDIRECTORY(Tutorial/ExitTest                  ${EXEC_WORKING_DIR}/Tutorial/ExitTest                 )
  ADD_SUBDIRECTORY(Tutorial/GeoCellMetrics_uns        ${EXEC_WORKING_DIR}/Tutorial/GeoCellMetrics_uns       )
  ADD_SUBDIRECTORY(Tutorial/GeoCellMetrics_str        ${EXEC_WORKING_DIR}/Tutorial/GeoCellMetrics_str       )
  ADD_SUBDIRECTORY(Tutorial/GeoFaceMetrics_uns        ${EXEC_WORKING_DIR}/Tutorial/GeoFaceMetrics_uns       )
  ADD_SUBDIRECTORY(Tutorial/GeoFaceMetrics_str        ${EXEC_WORKING_DIR}/Tutorial/GeoFaceMetrics_str       )
  ADD_SUBDIRECTORY(Tutorial/GeoNodeTopo_uns           ${EXEC_WORKING_DIR}/Tutorial/GeoNodeTopo_uns          )
  ADD_SUBDIRECTORY(Tutorial/GeoNodeTopo_str           ${EXEC_WORKING_DIR}/Tutorial/GeoNodeTopo_str          )
  ADD_SUBDIRECTORY(Tutorial/GeoFaceTopo_uns           ${EXEC_WORKING_DIR}/Tutorial/GeoFaceTopo_uns          )
  ADD_SUBDIRECTORY(Tutorial/GeoCellTopo_uns           ${EXEC_WORKING_DIR}/Tutorial/GeoCellTopo_uns          )
  ADD_SUBDIRECTORY(Tutorial/TimeSpanTest              ${EXEC_WORKING_DIR}/Tutorial/TimeSpanTest             )
  ADD_SUBDIRECTORY(Tutorial/DYutilTest                ${EXEC_WORKING_DIR}/Tutorial/DYutilTest               )
  ADD_SUBDIRECTORY(Tutorial/TK_ParseTest              ${EXEC_WORKING_DIR}/Tutorial/TK_ParseTest             )
  ADD_SUBDIRECTORY(Tutorial/CGNSConversionToPlot3D    ${EXEC_WORKING_DIR}/Tutorial/CGNSConversionToPlot3D   )
  ADD_SUBDIRECTORY(Tutorial/GlobalDataBaseTest        ${EXEC_WORKING_DIR}/Tutorial/GlobalDataBaseTest       )
  ADD_SUBDIRECTORY(Tutorial/SafeDataTest              ${EXEC_WORKING_DIR}/Tutorial/SafeDataTest             )
  ADD_SUBDIRECTORY(Tutorial/BC_str                    ${EXEC_WORKING_DIR}/Tutorial/BC_str                   )
  ADD_SUBDIRECTORY(Tutorial/BC_uns                    ${EXEC_WORKING_DIR}/Tutorial/BC_uns                   )
  ADD_SUBDIRECTORY(Tutorial/Data_UtilTest             ${EXEC_WORKING_DIR}/Tutorial/Data_UtilTest            )
  ADD_SUBDIRECTORY(Tutorial/Data_FieldTest            ${EXEC_WORKING_DIR}/Tutorial/Data_FieldTest           )
  ADD_SUBDIRECTORY(Tutorial/Data_IndexDataTest        ${EXEC_WORKING_DIR}/Tutorial/Data_IndexDataTest       )
  ADD_SUBDIRECTORY(Tutorial/Data_ParamFieldSuiteTest  ${EXEC_WORKING_DIR}/Tutorial/Data_ParamFieldSuiteTest )
  ADD_SUBDIRECTORY(Tutorial/FluentToPHengLEI_uns      ${EXEC_WORKING_DIR}/Tutorial/FluentToPHengLEI_uns     )
  #ADD_SUBDIRECTORY(Tutorial/PreGridMirror            ${EXEC_WORKING_DIR}/Tutorial                          )
  #ADD_SUBDIRECTORY(Tutorial/KDTree                   ${EXEC_WORKING_DIR}/Tutorial/KDTree                   )
  ADD_SUBDIRECTORY(Tutorial/WalldistCompute           ${EXEC_WORKING_DIR}/Tutorial/WalldistCompute          )
  #ADD_SUBDIRECTORY(Tutorial/UT_GridMemberFunction    ${EXEC_WORKING_DIR}/Tutorial/UT_GridMemberFunction    )
  
  #ADD_SUBDIRECTORY(Tutorial/SolverDemo ${EXEC_WORKING_DIR}/Tutorial)
endif(BUILD_TUTORIAL)
#exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
message(STATUS "Finished build project ${PROJECT_NAME}!")


