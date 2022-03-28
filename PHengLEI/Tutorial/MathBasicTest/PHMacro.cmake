#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                             +
#       Hybrid Platform for Engineering and Research of Flows (HyperFlow)     +
#                        £¨C£©Zhang Laiping and He Xin                        +
#                    State Key Laboratory of Aerodynamics                     +
#                    Computational Aerodynamics Institute                     +
#            China Aerodynamics Research and Development Center               +
#                                 Since 2007                                  +
#                                                                             +
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#! @file      PHMacro.cmake
#! @brief     It defines the macro of USE_XX precondition.
#! @author    Bell.

#Mac function to define pre-definition.
macro(opt OPTION HELP VALUE)
  option(USE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item USE_${OPTION}\n${HELP} (default: ${VALUE})")
  if(${VALUE} STREQUAL "ON")
    message(STATUS "${OPTION} is turned on by default")
  else()
	message(STATUS "${OPTION} is turned off by default")
  endif()
endmacro(opt)

macro(opt_build OPTION HELP VALUE)
  option(BUILD_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item BUILD_${OPTION}\n${HELP} (default: ${VALUE})")
  if(${VALUE} STREQUAL "ON")
    message(STATUS "${OPTION} is turned on by default")
  else()
	message(STATUS "${OPTION} is turned off by default")
  endif()
endmacro(opt_build)
