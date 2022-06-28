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

if(BUILD_PBSolver)
  message(STATUS "BUILD_PBSolver is turned on!")
  add_definitions(-DUSE_PBSOLVER)
endif()

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

# set directory of PBSolver:
if(BUILD_PBSolver)
set(PBSOLVER_H_PATH PBSolver/include)
set(PBSOLVER_H_LIST
${PBSOLVER_H_PATH}/PBBoundary.h
${PBSOLVER_H_PATH}/PBBoundaryFlow.h
${PBSOLVER_H_PATH}/PBBoundaryFlowPressureOutlet.h
${PBSOLVER_H_PATH}/PBBoundaryFlowSolidSurface.h
${PBSOLVER_H_PATH}/PBBoundaryFlowVelocityInlet.h
${PBSOLVER_H_PATH}/PBBoundaryScalar.h
${PBSOLVER_H_PATH}/PBBoundaryScalarConvection.h
${PBSOLVER_H_PATH}/PBBoundaryScalarDefault.h
${PBSOLVER_H_PATH}/PBBoundaryScalarDiffusion.h
${PBSOLVER_H_PATH}/PBBoundaryEnergy.h
${PBSOLVER_H_PATH}/PBBoundaryEnergyDefault.h
${PBSOLVER_H_PATH}/PBBoundaryEnergyInlet.h
${PBSOLVER_H_PATH}/PBBoundaryEnergyOutlet.h
${PBSOLVER_H_PATH}/PBBoundaryEnergyWall.h
${PBSOLVER_H_PATH}/PBBoundaryTurb.h
${PBSOLVER_H_PATH}/PBBoundaryTurbDefault.h
${PBSOLVER_H_PATH}/PBBoundaryTurbEpsilonInlet.h
${PBSOLVER_H_PATH}/PBBoundaryTurbEpsilonPressureOutlet.h
${PBSOLVER_H_PATH}/PBBoundaryTurbEpsilonWall.h
${PBSOLVER_H_PATH}/PBBoundaryTurbKInlet.h
${PBSOLVER_H_PATH}/PBBoundaryTurbKPressureOutlet.h
${PBSOLVER_H_PATH}/PBBoundaryTurbKWall.h
${PBSOLVER_H_PATH}/PBBoundaryTurbSAWall.h
${PBSOLVER_H_PATH}/PBBoundaryTurbSAInlet.h
${PBSOLVER_H_PATH}/PBBoundaryTurbSAPressureOutlet.h
${PBSOLVER_H_PATH}/PBConvection.h
${PBSOLVER_H_PATH}/PBConvectionCDS.h
${PBSOLVER_H_PATH}/PBConvectionQUICK.h
${PBSOLVER_H_PATH}/PBConvectionUDS.h
${PBSOLVER_H_PATH}/PBConvectionSUDS.h
${PBSOLVER_H_PATH}/PBDiffusion.h
${PBSOLVER_H_PATH}/PBDiffusionNonorthogonal.h
${PBSOLVER_H_PATH}/PBGrad.h
${PBSOLVER_H_PATH}/PBGradGauss.h
${PBSOLVER_H_PATH}/PBGradLsq.h
${PBSOLVER_H_PATH}/PBIterSolv.h
${PBSOLVER_H_PATH}/PBIterSolvBiCGSTAB.h
${PBSOLVER_H_PATH}/PBIterSolvGMRES.h
${PBSOLVER_H_PATH}/PBIterSolvAMG.h
${PBSOLVER_H_PATH}/PBObject.h
${PBSOLVER_H_PATH}/PBOutput.h
${PBSOLVER_H_PATH}/PBOutputTecplot.h
${PBSOLVER_H_PATH}/PBPreconditioner.h
${PBSOLVER_H_PATH}/PBPreconditionerDILU.h
${PBSOLVER_H_PATH}/PBS.h
${PBSOLVER_H_PATH}/PBSolver.h
${PBSOLVER_H_PATH}/PBSolverFlow.h
${PBSOLVER_H_PATH}/PBSolverScalar.h
${PBSOLVER_H_PATH}/PBSolverEnergy.h
${PBSOLVER_H_PATH}/PBSolverTurbK.h
${PBSOLVER_H_PATH}/PBSolverKEpsilon.h
${PBSOLVER_H_PATH}/PBSolverTurbSA.h
${PBSOLVER_H_PATH}/PBSolverTurbEpsilon.h
${PBSOLVER_H_PATH}/PBSource.h
${PBSOLVER_H_PATH}/PBSourceFlowDefault.h
${PBSOLVER_H_PATH}/PBSourceFlowGravity.h
${PBSOLVER_H_PATH}/PBSourceEnergyDefault.h
${PBSOLVER_H_PATH}/PBSourceEnergyDissipation.h
${PBSOLVER_H_PATH}/PBSourceTurbEpsilonDefault.h
${PBSOLVER_H_PATH}/PBSourceTurbKDefault.h
${PBSOLVER_H_PATH}/PBSourceTurbSADefault.h
${PBSOLVER_H_PATH}/PBTransient.h
${PBSOLVER_H_PATH}/PBTransientCrankNicolson.h
${PBSOLVER_H_PATH}/PBTransientImplicit2ndOrder.h
${PBSOLVER_H_PATH}/PBTransientImplicitEuler.h
)
source_group(PBSolver\\include FILES ${PBSOLVER_H_LIST})
include_directories(${PBSOLVER_H_PATH})

set(PBSOLVER_CPP_PATH PBSolver/src)
set(PBSOLVER_CPP_LIST
${PBSOLVER_CPP_PATH}/PBBoundary.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryFlow.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryFlowPressureOutlet.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryFlowSolidSurface.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryFlowVelocityInlet.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryScalar.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryScalarConvection.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryScalarDefault.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryScalarDiffusion.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryEnergy.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryEnergyDefault.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryEnergyInlet.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryEnergyOutlet.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryEnergyWall.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryTurb.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryTurbDefault.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryTurbEpsilonInlet.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryTurbEpsilonPressureOutlet.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryTurbEpsilonWall.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryTurbSAWall.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryTurbSAInlet.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryTurbSAPressureOutlet.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryTurbKInlet.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryTurbKPressureOutlet.cpp
${PBSOLVER_CPP_PATH}/PBBoundaryTurbKWall.cpp
${PBSOLVER_CPP_PATH}/PBConvection.cpp
${PBSOLVER_CPP_PATH}/PBConvectionCDS.cpp
${PBSOLVER_CPP_PATH}/PBConvectionQUICK.cpp
${PBSOLVER_CPP_PATH}/PBConvectionUDS.cpp
${PBSOLVER_CPP_PATH}/PBConvectionSUDS.cpp
${PBSOLVER_CPP_PATH}/PBDiffusion.cpp
${PBSOLVER_CPP_PATH}/PBDiffusionNonorthogonal.cpp
${PBSOLVER_CPP_PATH}/PBGrad.cpp
${PBSOLVER_CPP_PATH}/PBGradGauss.cpp
${PBSOLVER_CPP_PATH}/PBGradLsq.cpp
${PBSOLVER_CPP_PATH}/PBIterSolv.cpp
${PBSOLVER_CPP_PATH}/PBIterSolvBiCGSTAB.cpp
${PBSOLVER_CPP_PATH}/PBIterSolvGMRES.cpp
${PBSOLVER_CPP_PATH}/PBIterSolvAMG.cpp
${PBSOLVER_CPP_PATH}/PBOutput.cpp
${PBSOLVER_CPP_PATH}/PBOutputTecplot.cpp
${PBSOLVER_CPP_PATH}/PBPreconditioner.cpp
${PBSOLVER_CPP_PATH}/PBPreconditionerDILU.cpp
${PBSOLVER_CPP_PATH}/PBSolver.cpp
${PBSOLVER_CPP_PATH}/PBSolverFlow.cpp
${PBSOLVER_CPP_PATH}/PBSolverScalar.cpp
${PBSOLVER_CPP_PATH}/PBSolverEnergy.cpp
${PBSOLVER_CPP_PATH}/PBSolverTurbEpsilon.cpp
${PBSOLVER_CPP_PATH}/PBSolverTurbK.cpp
${PBSOLVER_CPP_PATH}/PBSolverKEpsilon.cpp
${PBSOLVER_CPP_PATH}/PBSolverTurbSA.cpp
${PBSOLVER_CPP_PATH}/PBSource.cpp
${PBSOLVER_CPP_PATH}/PBSourceFlowDefault.cpp
${PBSOLVER_CPP_PATH}/PBSourceFlowGravity.cpp
${PBSOLVER_CPP_PATH}/PBSourceEnergyDefault.cpp
${PBSOLVER_CPP_PATH}/PBSourceEnergyDissipation.cpp
${PBSOLVER_CPP_PATH}/PBSourceTurbEpsilonDefault.cpp
${PBSOLVER_CPP_PATH}/PBSourceTurbKDefault.cpp
${PBSOLVER_CPP_PATH}/PBSourceTurbSADefault.cpp
${PBSOLVER_CPP_PATH}/PBTransient.cpp
${PBSOLVER_CPP_PATH}/PBTransientCrankNicolson.cpp
${PBSOLVER_CPP_PATH}/PBTransientImplicit2ndOrder.cpp
${PBSOLVER_CPP_PATH}/PBTransientImplicitEuler.cpp
)
source_group(PBSolver\\src FILES ${PBSOLVER_CPP_LIST})
endif(BUILD_PBSolver)

#20.1 set header file list
set(CFD_HEAD_LIST 
${SOLVER_H_LIST}
${NSSOLVER_H_LIST}
${TurbSOLVER_H_LIST}
${SpecDiffHybSolver_H_LIST}
${HOSOLVERUNSTRUCT_H_LIST}
${DEMOSOLVER_H_LIST}
${PBSOLVER_H_LIST}
)

# add INHOUSE_VERSION_H_LIST

# add SpecDiffHybSolver_H_LIST
if(BUILD_SpecDiffHybSolver)
list(APPEND CFD_HEAD_LIST ${SpecDiffHybSolver_H_LIST})
endif()

# add DEMOSOLVER_H_LIST
if(BUILD_DEMOSOLVER)
list(APPEND CFD_HEAD_LIST ${DEMOSOLVER_H_LIST})
endif()

# add PBSOLVER_H_LIST
if(BUILD_PBSolver)
list(APPEND CFD_HEAD_LIST ${PBSOLVER_H_LIST})
endif()

#20.2 set source file list:
set(CFD_CPP_LIST 
${SOLVER_CPP_LIST}
${NSSOLVER_CPP_LIST}
${TurbSOLVER_CPP_LIST}
${SpecDiffHybSolver_CPP_LIST}
${HOSOLVERUNSTRUCT_CPP_LIST}
${DEMOSOLVER_CPP_LIST}
${PBSOLVER_CPP_LIST}
)
# add INHOUSE_VERSION_CPP_LIST

# add SpecDiffHybSolver_CPP_LIST
if(BUILD_SpecDiffHybSolver)
list(APPEND CFD_CPP_LIST ${SpecDiffHybSolver_CPP_LIST})
endif()

# add DEMOSOLVER_CPP_LIST
if(BUILD_DEMOSOLVER)
list(APPEND CFD_CPP_LIST ${DEMOSOLVER_CPP_LIST})
endif()

# add PBSOLVER_CPP_LIST
if(BUILD_PBSolver)
list(APPEND CFD_CPP_LIST ${PBSOLVER_CPP_LIST})
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

if(WIN32)
  list(APPEND LINK_LIBRARY_LIST msmpi)
endif(WIN32)

if(BUILD_SpecDiffHybSolver)
  list(APPEND LINK_LIBRARY_LIST p3dfft fftw3 mpifort gfortran)
endif(BUILD_SpecDiffHybSolver)

target_link_libraries(${CFD_TARGET_NAME} ${LINK_LIBRARY_LIST})

message(STATUS "Link library list: ${LINK_LIBRARY_LIST}.")


# For OpenMP
if(USE_OMP)
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
endif(USE_OMP)

#24. exit:
message(STATUS "--- Exit ${CMAKE_CURRENT_LIST_FILE}! ---")

