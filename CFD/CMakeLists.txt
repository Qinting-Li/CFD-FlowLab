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
include(${PROJECT_SOURCE_DIR}/../PHengLEI/PHMacro.cmake)
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
if(USE_INHOUSE_VERSION)
  message(STATUS "USE_INHOUSE_VERSION is turned on!")
  add_definitions(-DUSE_INHOUSE_VERSION)
endif(USE_INHOUSE_VERSION)
if(WIN32)
if(USE_Windows_X64)
  message(STATUS "USE_Windows_X64 is turned on!")
  add_definitions(-DUSE_WINDOWS_X64)
endif(USE_Windows_X64)
endif(WIN32)
if(BUILD_SpecDiffHybSolver)
  message(STATUS "BUILD_SpecDiffHybSolver is turned on!")
  add_definitions(-DUSE_SpecDiffHybSolver)
endif()
if(BUILD_IncompressibleSolver)
  message(STATUS "BUILD_IncompressibleSolver is turned on!")
  add_definitions(-DUSE_IncompressibleSolver)
endif(BUILD_IncompressibleSolver)
#1. output message:
message(STATUS "") 
message(STATUS "--- Enter subdirectory <CFD>, Start to call ${CMAKE_CURRENT_LIST_FILE} ---") 
message(STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}.")
message(STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}.")

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
link_directories("${LIB_API_EXPORT_PATH}")

#if(WIN32)
#if(USE_Windows_X64)
#  set(LIB_HYFLOW_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../HyFLOW/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
#  set(CMAKE_EXE_LINKER_FLAGS "/machine:x64")
#else(USE_Windows_X64)
#  set(LIB_HYFLOW_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../HyFLOW/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
#endif(USE_Windows_X64)
#else(WIN32)
#  set(LIB_HYFLOW_EXPORT_PATH ${PROJECT_SOURCE_DIR}/../HyFLOW/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
#endif(WIN32)
#link_directories("${LIB_HYFLOW_EXPORT_PATH}")

set(CFD_TARGET_NAME "PHCFD")
set(LIBRARY_OUTPUT_PATH)
if(WIN32)
if(USE_Windows_X64)
 set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE})
 set(CMAKE_EXE_LINKER_FLAGS "/machine:x64")
else(USE_Windows_X64)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE})
endif(USE_Windows_X64)
else(WIN32)
  set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../CFD/bin/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
endif(WIN32)
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

includeHeadFiles()

set(HDF5_H_PATH    ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/HDF5/include)
include_directories("${HDF5_H_PATH}")

if(USE_TecplotLib)
set(TECIO_H_PATH   ${PROJECT_SOURCE_DIR}/../PHengLEI/3rdparty/TECIO/include)
include_directories("${TECIO_H_PATH}")
endif(USE_TecplotLib)

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
${NSSOLVER_H_PATH}/MultiGridOperation.h
${NSSOLVER_H_PATH}/Param_NSSolver.h
${NSSOLVER_H_PATH}/Param_NSSolverUnstruct.h
${NSSOLVER_H_PATH}/Param_NSSolverStruct.h
${NSSOLVER_H_PATH}/HyFlowSolverStruct.h
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
${NSSOLVER_CPP_PATH}/Param_NSSolver.cpp
${NSSOLVER_CPP_PATH}/Param_NSSolverUnstruct.cpp
${NSSOLVER_CPP_PATH}/Param_NSSolverStruct.cpp
${NSSOLVER_CPP_PATH}/HyFlowSolverStruct.cpp
)
source_group(NSSolver\\src FILES ${NSSOLVER_CPP_LIST})

# set directory of TurbSolver:
set(TurbSOLVER_H_PATH        TurbSolver/include)
set(TurbSOLVER_H_LIST
${TurbSOLVER_H_PATH}/TurbSolver.h
${TurbSOLVER_H_PATH}/TurbSolverStruct.h
${TurbSOLVER_H_PATH}/TurbSolverUnstr.h
${TurbSOLVER_H_PATH}/Transition.h
${TurbSOLVER_H_PATH}/LESSolver.h
${TurbSOLVER_H_PATH}/LESSolverStruct.h
${TurbSOLVER_H_PATH}/LESSolverUnstr.h
${TurbSOLVER_H_PATH}/TurbSolverStructFD.h
${TurbSOLVER_H_PATH}/Param_TurbSolver.h
${TurbSOLVER_H_PATH}/Param_TurbSolverStruct.h
${TurbSOLVER_H_PATH}/Param_TurbSolverUnstruct.h
${TurbSOLVER_H_PATH}/Param_LESSolver.h
${TurbSOLVER_H_PATH}/Param_LESSolverStruct.h
${TurbSOLVER_H_PATH}/Param_LESSolverUnstruct.h
)
source_group(TurbSolver\\include FILES ${TurbSOLVER_H_LIST}) 
include_directories(${TurbSOLVER_H_PATH})

set(TurbSOLVER_CPP_PATH      TurbSolver/src)
set(TurbSOLVER_CPP_LIST
${TurbSOLVER_CPP_PATH}/TurbSolver.cpp
${TurbSOLVER_CPP_PATH}/TurbSolverStruct.cpp
${TurbSOLVER_CPP_PATH}/TurbSolverStructFD.cpp
${TurbSOLVER_CPP_PATH}/TurbSolverUnstr.cpp
${TurbSOLVER_CPP_PATH}/Transition.cpp
${TurbSOLVER_CPP_PATH}/LESSolver.cpp
${TurbSOLVER_CPP_PATH}/LESSolverStruct.cpp
${TurbSOLVER_CPP_PATH}/LESSolverUnstr.cpp
${TurbSOLVER_CPP_PATH}/Param_TurbSolver.cpp
${TurbSOLVER_CPP_PATH}/Param_TurbSolverStruct.cpp
${TurbSOLVER_CPP_PATH}/Param_TurbSolverUnstruct.cpp
${TurbSOLVER_CPP_PATH}/Param_LESSolver.cpp
${TurbSOLVER_CPP_PATH}/Param_LESSolverStruct.cpp
${TurbSOLVER_CPP_PATH}/Param_LESSolverUnstruct.cpp
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

if(BUILD_LBM)
set(LBMSOLVER_H_PATH           LBMSolver/include)
set(LBMSOLVER_H_LIST
${LBMSOLVER_H_PATH}/BinBase.h
${LBMSOLVER_H_PATH}/BoundaryCondition.h
${LBMSOLVER_H_PATH}/BoundaryFactory.h
${LBMSOLVER_H_PATH}/cavity.h
${LBMSOLVER_H_PATH}/Collision.h
${LBMSOLVER_H_PATH}/CommonAPI.h
${LBMSOLVER_H_PATH}/CommonLibDll.h
${LBMSOLVER_H_PATH}/D2Q9.h
${LBMSOLVER_H_PATH}/D3Q15.h
${LBMSOLVER_H_PATH}/D3Q19.h
${LBMSOLVER_H_PATH}/DevEnv.h
${LBMSOLVER_H_PATH}/DimClass.h
${LBMSOLVER_H_PATH}/Function.h
${LBMSOLVER_H_PATH}/GenBcMap.h
${LBMSOLVER_H_PATH}/GenColorMap.h
${LBMSOLVER_H_PATH}/Geo3d.h
${LBMSOLVER_H_PATH}/GetConfig.h
${LBMSOLVER_H_PATH}/GlobalData.h
${LBMSOLVER_H_PATH}/InitField.h
${LBMSOLVER_H_PATH}/json.hpp
${LBMSOLVER_H_PATH}/JsonGeoGen.h
${LBMSOLVER_H_PATH}/JsonGridGen.h
${LBMSOLVER_H_PATH}/JsonRangeGen.h
${LBMSOLVER_H_PATH}/JsonReader.h
${LBMSOLVER_H_PATH}/Lattice.h
${LBMSOLVER_H_PATH}/LatticeMPI.h
${LBMSOLVER_H_PATH}/LBMGeometry.h
${LBMSOLVER_H_PATH}/MemPool.h
${LBMSOLVER_H_PATH}/PHLBM.h
${LBMSOLVER_H_PATH}/PHLBMGenerator.h
${LBMSOLVER_H_PATH}/PHLBMGeoGen.h
${LBMSOLVER_H_PATH}/PHLBMGridGen.h
${LBMSOLVER_H_PATH}/PHLBMRangeGen.h
${LBMSOLVER_H_PATH}/PHLBMUnstructGrid.h
${LBMSOLVER_H_PATH}/Range.h
${LBMSOLVER_H_PATH}/StlGeo.h
${LBMSOLVER_H_PATH}/TimeMarching.h
${LBMSOLVER_H_PATH}/UniformGrid.h
)
source_group(LBMSolver\\include FILES ${LBMSOLVER_H_LIST}) 
include_directories(${LBMSOLVER_H_PATH})

set(LBMSOLVER_CPP_PATH      LBMSolver/src)
set(LBMSOLVER_CPP_LIST
${LBMSOLVER_CPP_PATH}/BinBase.cpp
${LBMSOLVER_CPP_PATH}/cavity-mpi.cpp
${LBMSOLVER_CPP_PATH}/Collision.cpp
${LBMSOLVER_CPP_PATH}/D2Q9.cpp
${LBMSOLVER_CPP_PATH}/D3Q15.cpp
${LBMSOLVER_CPP_PATH}/D3Q19.cpp
${LBMSOLVER_CPP_PATH}/EnrollBinReader.cpp
${LBMSOLVER_CPP_PATH}/Function.cpp
${LBMSOLVER_CPP_PATH}/Geometry.cpp
${LBMSOLVER_CPP_PATH}/GetConfig.cpp
${LBMSOLVER_CPP_PATH}/GlobalData.cpp
${LBMSOLVER_CPP_PATH}/JsonReader.cpp
${LBMSOLVER_CPP_PATH}/LatticeMPI.cpp
${LBMSOLVER_CPP_PATH}/PHLBMBCGen.cpp
${LBMSOLVER_CPP_PATH}/PHLBMGeoGen.cpp
${LBMSOLVER_CPP_PATH}/PHLBMUnstructGrid.cpp
${LBMSOLVER_CPP_PATH}/StlGeo.cpp
)
source_group(LBMSolver\\src FILES ${LBMSOLVER_CPP_LIST})
endif(BUILD_LBM)
# Build In-house version:
if(BUILD_IncompressibleSolver)
# set directory of IncompressibleSolver:
set(IncompressibleSolver_H_PATH        IncompressibleSolver/include)
set(IncompressibleSolver_H_LIST
${IncompressibleSolver_H_PATH}/INCompSolver.h
${IncompressibleSolver_H_PATH}/INCompSolverStruct.h
${IncompressibleSolver_H_PATH}/INCompSolverUnstruct.h
${IncompressibleSolver_H_PATH}/Param_INCompSolver.h
${IncompressibleSolver_H_PATH}/Param_INCompSolverStruct.h
${IncompressibleSolver_H_PATH}/Param_INCompSolverUnstruct.h
${IncompressibleSolver_H_PATH}/DiffCoeffCalculator.h
${IncompressibleSolver_H_PATH}/DiffusionCoeff.h
${IncompressibleSolver_H_PATH}/EquationOfState.h
${IncompressibleSolver_H_PATH}/InCom_Calculator.h
${IncompressibleSolver_H_PATH}/InCom_Energy_Calculator.h
${IncompressibleSolver_H_PATH}/InCom_Scalar_Calculator.h
${IncompressibleSolver_H_PATH}/InCom_Turbulence_Calculator.h
${IncompressibleSolver_H_PATH}/InComEffectiveViscCalculator.h
${IncompressibleSolver_H_PATH}/InComEnergyDiffusionCoefficient.h
${IncompressibleSolver_H_PATH}/VFluxCalculator.h
${IncompressibleSolver_H_PATH}/ViscLamCalculator.h
${IncompressibleSolver_H_PATH}/PPESIMPLESolver.h
${IncompressibleSolver_H_PATH}/PPEProjectionSolver.h
${IncompressibleSolver_H_PATH}/ScalarSolver.h
${IncompressibleSolver_H_PATH}/PPESolver.h
${IncompressibleSolver_H_PATH}/LinearEquationSolver.h
${IncompressibleSolver_H_PATH}/MomEquationSolver.h
)
source_group(IncompressibleSolver\\include FILES ${IncompressibleSolver_H_LIST}) 
include_directories(${IncompressibleSolver_H_PATH})

set(IncompressibleSolver_CPP_PATH      IncompressibleSolver/src)
set(IncompressibleSolver_CPP_LIST
${IncompressibleSolver_CPP_PATH}/INCompSolver.cpp
${IncompressibleSolver_CPP_PATH}/INCompSolverStruct.cpp
${IncompressibleSolver_CPP_PATH}/INCompSolverUnstruct.cpp
${IncompressibleSolver_CPP_PATH}/Param_INCompSolver.cpp
${IncompressibleSolver_CPP_PATH}/Param_INCompSolverStruct.cpp
${IncompressibleSolver_CPP_PATH}/Param_INCompSolverUnstruct.cpp
${IncompressibleSolver_CPP_PATH}/DiffCoeffCalculator.cpp
${IncompressibleSolver_CPP_PATH}/DiffusionCoeff.cpp
${IncompressibleSolver_CPP_PATH}/EquationOfState.cpp
${IncompressibleSolver_CPP_PATH}/InCom_Calculator.cpp
${IncompressibleSolver_CPP_PATH}/InCom_Energy_Calculator.cpp
${IncompressibleSolver_CPP_PATH}/InCom_Scalar_Calculator.cpp
${IncompressibleSolver_CPP_PATH}/InCom_Turbulence_Calculator.cpp
${IncompressibleSolver_CPP_PATH}/InComEffectiveViscCalculator.cpp
${IncompressibleSolver_CPP_PATH}/InComEnergyDiffusionCoefficient.cpp
${IncompressibleSolver_CPP_PATH}/VFluxCalculator.cpp
${IncompressibleSolver_CPP_PATH}/ViscLamCalculator.cpp
${IncompressibleSolver_CPP_PATH}/PPESIMPLESolver.cpp
${IncompressibleSolver_CPP_PATH}/PPEProjectionSolver.cpp
${IncompressibleSolver_CPP_PATH}/MomEquationSolver.cpp
${IncompressibleSolver_CPP_PATH}/ScalarSolver.cpp
${IncompressibleSolver_CPP_PATH}/PPESolver.cpp
${IncompressibleSolver_CPP_PATH}/LinearEquationSolver.cpp
)
source_group(IncompressibleSolver\\src FILES ${IncompressibleSolver_CPP_LIST})

endif(BUILD_IncompressibleSolver)

if(USE_INHOUSE_VERSION)
# set directory of DSMCSolver:
set(DSMCSOLVER_H_PATH        DSMCSolver/include)
set(DSMCSOLVER_H_LIST
${DSMCSOLVER_H_PATH}/DSMCSolver.h
${DSMCSOLVER_H_PATH}/DSMCCollision.h
${DSMCSOLVER_H_PATH}/DSMCPointers.h
${DSMCSOLVER_H_PATH}/DSMCMolecule.h
${DSMCSOLVER_H_PATH}/DSMCGrid.h
${DSMCSOLVER_H_PATH}/DSMCConsts.h
${DSMCSOLVER_H_PATH}/DSMCSample.h
${DSMCSOLVER_H_PATH}/DSMCIndex.h
)
source_group(DSMCSolver\\include FILES ${DSMCSOLVER_H_LIST}) 
include_directories(${DSMCSOLVER_H_PATH})

set(DSMCSOLVER_CPP_PATH      DSMCSolver/src)
set(DSMCSOLVER_CPP_LIST
${DSMCSOLVER_CPP_PATH}/DSMCSolver.cpp
${DSMCSOLVER_CPP_PATH}/DSMCCollision.cpp
${DSMCSOLVER_CPP_PATH}/DSMCMolecule.cpp
${DSMCSOLVER_CPP_PATH}/DSMCGrid.cpp
${DSMCSOLVER_CPP_PATH}/DSMCSample.cpp
${DSMCSOLVER_CPP_PATH}/DSMCIndex.cpp
)
source_group(DSMCSolver\\src FILES ${DSMCSOLVER_CPP_LIST})

# set directory of Helicopter module
set(HelicopterModule_H_PATH HelicopterModule/include)
set(HelicopterModule_H_LIST
${HelicopterModule_H_PATH}/MBConstTable.h
${HelicopterModule_H_PATH}/HelicopterSolver.h
${HelicopterModule_H_PATH}/RotaryCoordinateSystem.h
${HelicopterModule_H_PATH}/MBGridFactory.h
${HelicopterModule_H_PATH}/MBGridGroup.h
${HelicopterModule_H_PATH}/MBModualCenter.h
${HelicopterModule_H_PATH}/MBReactor.h
)
source_group(HelicopterModule\\include FILES ${HelicopterModule_H_LIST})
include_directories(${HelicopterModule_H_PATH})

set(HelicopterModule_CPP_PATH HelicopterModule/src)
set(HelicopterModule_CPP_LIST
${HelicopterModule_CPP_PATH}/HelicopterSolver.cpp
${HelicopterModule_CPP_PATH}/RotaryCoordinateSystem.cpp
${HelicopterModule_CPP_PATH}/MBGridFactory.cpp
${HelicopterModule_CPP_PATH}/MBGridGroup.cpp
${HelicopterModule_CPP_PATH}/MBModualCenter.cpp
${HelicopterModule_CPP_PATH}/MBReactor.cpp
)
source_group(HelicopterModule\\src FILES ${HelicopterModule_CPP_LIST})

# set directory of UGKSSolver:
set(UGKSSOLVER_H_PATH        UGKSSolver/include)
set(UGKSSOLVER_H_LIST
${UGKSSOLVER_H_PATH}/UGKSSolver.h
${UGKSSOLVER_H_PATH}/UGKSSimulation.h
${UGKSSOLVER_H_PATH}/Param_UGKSSolver.h
${UGKSSOLVER_H_PATH}/UGKSContants.h
${UGKSSOLVER_H_PATH}/UGKSLimiter.h
${UGKSSOLVER_H_PATH}/FluxParamProxy.h
)
source_group(UGKSSolver\\include FILES ${UGKSSOLVER_H_LIST}) 
include_directories(${UGKSSOLVER_H_PATH})

set(UGKSSOLVER_CPP_PATH      UGKSSolver/src)
set(UGKSSOLVER_CPP_LIST
${UGKSSOLVER_CPP_PATH}/UGKSSolver.cpp
${UGKSSOLVER_CPP_PATH}/UGKSSimulation.cpp
${UGKSSOLVER_CPP_PATH}/Param_UGKSSolver.cpp
${UGKSSOLVER_CPP_PATH}/FluxParamProxy.cpp
)
source_group(UGKSSolver\\src FILES ${UGKSSOLVER_CPP_LIST})
endif(USE_INHOUSE_VERSION)

if(USE_INHOUSE_VERSION)
if(WIN32)
if(USE_Windows_X64)
    LINK_DIRECTORIES("${PROJECT_SOURCE_DIR}/../HelicopterModule/lib/${CMAKE_SYSTEM_NAME}-x64/${CMAKE_BUILD_TYPE}")
else(USE_Windows_X64)
    LINK_DIRECTORIES("${PROJECT_SOURCE_DIR}/../HelicopterModule/lib/${CMAKE_SYSTEM_NAME}-x86/${CMAKE_BUILD_TYPE}")
endif(USE_Windows_X64)
else(WIN32)
    LINK_DIRECTORIES("${PROJECT_SOURCE_DIR}/../HelicopterModule/lib/Linux")
endif(WIN32)
endif(USE_INHOUSE_VERSION)

# set directory of SpecDiffHybSolver:
if(BUILD_SpecDiffHybSolver)
set(SpecDiffHybSolver_H_PATH SpecDiffHybSolver/include)
set(SpecDiffHybSolver_H_LIST
${SpecDiffHybSolver_H_PATH}/SpecDiffHybSolver.h
${SpecDiffHybSolver_H_PATH}/SpecGrid.h
${SpecDiffHybSolver_H_PATH}/SpecDiffHybGrid.h
${SpecDiffHybSolver_H_PATH}/ExplicitDifferenceBoundary.h
${SpecDiffHybSolver_H_PATH}/CompactDifferenceFirstDerivative.h
${SpecDiffHybSolver_H_PATH}/CompactDifferenceSecondDerivative.h
${SpecDiffHybSolver_H_PATH}/PoissonSolver.h
${SpecDiffHybSolver_H_PATH}/CorrectWallDivergence.h
${SpecDiffHybSolver_H_PATH}/Complex.h
${SpecDiffHybSolver_H_PATH}/RescaleField.h
${SpecDiffHybSolver_H_PATH}/Statistics.h
${SpecDiffHybSolver_H_PATH}/p3dfft.h
${SpecDiffHybSolver_H_PATH}/fftw3.h
${SpecDiffHybSolver_H_PATH}/SpecDiffHybFieldView.h
${SpecDiffHybSolver_H_PATH}/LES.h
${SpecDiffHybSolver_H_PATH}/SpecSolver.h
${SpecDiffHybSolver_H_PATH}/Spectrum.h
${SpecDiffHybSolver_H_PATH}/Param_SpecSolver.h
${SpecDiffHybSolver_H_PATH}/Param_SpecDiffHybSolver.h
)
source_group(SpecDiffHybSolver\\include FILES ${SpecDiffHybSolver_H_LIST})
include_directories(${SpecDiffHybSolver_H_PATH})

set(SpecDiffHybSolver_CPP_PATH SpecDiffHybSolver/src)
set(SpecDiffHybSolver_CPP_LIST
${SpecDiffHybSolver_CPP_PATH}/SpecDiffHybSolver.cpp
${SpecDiffHybSolver_CPP_PATH}/SpecGrid.cpp
${SpecDiffHybSolver_CPP_PATH}/SpecDiffHybGrid.cpp
${SpecDiffHybSolver_CPP_PATH}/ExplicitDifferenceBoundary.cpp
${SpecDiffHybSolver_CPP_PATH}/CompactDifferenceFirstDerivative.cpp
${SpecDiffHybSolver_CPP_PATH}/CompactDifferenceSecondDerivative.cpp
${SpecDiffHybSolver_CPP_PATH}/PoissonSolver.cpp
${SpecDiffHybSolver_CPP_PATH}/CorrectWallDivergence.cpp
${SpecDiffHybSolver_CPP_PATH}/RescaleField.cpp
${SpecDiffHybSolver_CPP_PATH}/Statistics.cpp
${SpecDiffHybSolver_CPP_PATH}/SpecDiffHybFieldView.cpp
${SpecDiffHybSolver_CPP_PATH}/LES.cpp
${SpecDiffHybSolver_CPP_PATH}/SpecSolver.cpp
${SpecDiffHybSolver_CPP_PATH}/Spectrum.cpp
${SpecDiffHybSolver_CPP_PATH}/Param_SpecSolver.cpp
${SpecDiffHybSolver_CPP_PATH}/Param_SpecDiffHybSolver.cpp
)
source_group(SpecDiffHybSolver\\src FILES ${SpecDiffHybSolver_CPP_LIST})
endif(BUILD_SpecDiffHybSolver)

# set directory of DemoSolver:
if(BUILD_DEMOSOLVER)
set(DEMOSOLVER_H_PATH DemoSolver/include)
set(DEMOSOLVER_H_LIST
${DEMOSOLVER_H_PATH}/DemoSolver.h
${DEMOSOLVER_H_PATH}/Param_DemoSolver.h
)
source_group(DemoSolver\\include FILES ${DEMOSOLVER_H_LIST})
include_directories(${DEMOSOLVER_H_PATH})

set(DEMOSOLVER_CPP_PATH DemoSolver/src)
set(DEMOSOLVER_CPP_LIST
${DEMOSOLVER_CPP_PATH}/DemoSolver.cpp
${DEMOSOLVER_CPP_PATH}/Param_DemoSolver.cpp
)
source_group(DemoSolver\\src FILES ${DEMOSOLVER_CPP_LIST})
endif(BUILD_DEMOSOLVER)

#20.1 set header file list
set(CFD_HEAD_LIST 
${SOLVER_H_LIST}
${NSSOLVER_H_LIST}
${TurbSOLVER_H_LIST}
${SpecDiffHybSolver_H_LIST}
${HOSOLVERUNSTRUCT_H_LIST}
${LBMSOLVER_H_LIST}
${DEMOSOLVER_H_LIST}
)

# add INHOUSE_VERSION_H_LIST
if(USE_INHOUSE_VERSION)
list(APPEND CFD_HEAD_LIST ${HelicopterModule_H_LIST} ${DSMCSOLVER_H_LIST} ${UGKSSOLVER_H_LIST})
endif(USE_INHOUSE_VERSION)

# add SpecDiffHybSolver_H_LIST
if(BUILD_SpecDiffHybSolver)
list(APPEND CFD_HEAD_LIST ${SpecDiffHybSolver_H_LIST})
endif()

# add IncompressibleSolver_H_LIST
if(BUILD_IncompressibleSolver)
list(APPEND CFD_HEAD_LIST ${IncompressibleSolver_H_LIST})
endif(BUILD_IncompressibleSolver)

# add DEMOSOLVER_H_LIST
if(BUILD_DEMOSOLVER)
list(APPEND CFD_HEAD_LIST ${DEMOSOLVER_H_LIST})
endif()

#20.2 set source file list:
set(CFD_CPP_LIST 
${SOLVER_CPP_LIST}
${NSSOLVER_CPP_LIST}
${TurbSOLVER_CPP_LIST}
${SpecDiffHybSolver_CPP_LIST}
${HOSOLVERUNSTRUCT_CPP_LIST}
${LBMSOLVER_CPP_LIST}
${DEMOSOLVER_CPP_LIST}

)
# add INHOUSE_VERSION_CPP_LIST
if(USE_INHOUSE_VERSION)
list(APPEND CFD_CPP_LIST ${HelicopterModule_CPP_LIST} ${DSMCSOLVER_CPP_LIST} ${UGKSSOLVER_CPP_LIST})
endif(USE_INHOUSE_VERSION)

# add SpecDiffHybSolver_CPP_LIST
if(BUILD_SpecDiffHybSolver)
list(APPEND CFD_CPP_LIST ${SpecDiffHybSolver_CPP_LIST})
endif()

# add IncompressibleSolver_CPP_LIST
if(BUILD_IncompressibleSolver)
list(APPEND CFD_CPP_LIST ${IncompressibleSolver_CPP_LIST})
endif(BUILD_IncompressibleSolver)

# add DEMOSOLVER_CPP_LIST
if(BUILD_DEMOSOLVER)
list(APPEND CFD_CPP_LIST ${DEMOSOLVER_CPP_LIST})
endif()
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
set(API_TARGET_NAME "PHAPI")
list(APPEND LINK_LIBRARY_LIST ${API_TARGET_NAME} metis parmetis cgns hdf5)

if(USE_TecplotLib)
  list(APPEND LINK_LIBRARY_LIST tecio)
endif(USE_TecplotLib)

if(WIN32)
  list(APPEND LINK_LIBRARY_LIST msmpi)
endif(WIN32)

if(BUILD_SpecDiffHybSolver)
  list(APPEND LINK_LIBRARY_LIST p3dfft fftw3 mpifort gfortran)
endif(BUILD_SpecDiffHybSolver)

target_link_libraries(${CFD_TARGET_NAME} ${LINK_LIBRARY_LIST})

message(STATUS "Link library list: ${LINK_LIBRARY_LIST}.")


# For OpenMP
find_package (OpenMP)
if (OpenMP_FOUND)

if(WIN32)
target_include_directories (${CFD_TARGET_NAME}
    PRIVATE
        ${OpenMP_CXX_INCLUDE_DIRS}
)

target_compile_options (${CFD_TARGET_NAME}
    PRIVATE
        ${OpenMP_CXX_FLAGS}
)

target_compile_features (${CFD_TARGET_NAME}
    PUBLIC 
        cxx_std_14
)

target_link_libraries (${CFD_TARGET_NAME}
    PRIVATE 
        ${OpenMP_CXX_LIBRARIES}
)
else(WIN32) 
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(WIN32)
endif (OpenMP_FOUND)

#24. exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")
