# Copyright 1998-2019 Lawrence Livermore National Security, LLC and other
# HYPRE Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was HYPREConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include(CMakeFindDependencyMacro)

set(HYPRE_ENABLE_SHARED OFF)
set(HYPRE_ENABLE_BIGINT OFF)
set(HYPRE_ENABLE_MIXEDINT OFF)
set(HYPRE_ENABLE_SINGLE OFF)
set(HYPRE_ENABLE_LONG_DOUBLE OFF)
set(HYPRE_ENABLE_COMPLEX OFF)
set(HYPRE_ENABLE_HYPRE_BLAS ON)
set(HYPRE_ENABLE_HYPRE_LAPACK ON)
set(HYPRE_ENABLE_PERSISTENT_COMM OFF)
set(HYPRE_ENABLE_FEI OFF)
set(HYPRE_WITH_MPI ON)
set(HYPRE_WITH_OPENMP OFF)
set(HYPRE_WITH_HOPSCOTCH )
set(HYPRE_USING_DSUPERLU )
set(HYPRE_WITH_CALIPER OFF)
set(HYPRE_PRINT_ERRORS OFF)
set(HYPRE_TIMING OFF)
set(HYPRE_BUILD_EXAMPLES OFF)
set(HYPRE_BUILD_TESTS OFF)
set(HYPRE_USING_HOST_MEMORY ON)
set(HYPRE_WITH_CUDA OFF)
set(HYPRE_ENABLE_UNIFIED_MEMORY OFF)
set(HYPRE_ENABLE_CUDA_STREAMS ON)
set(HYPRE_ENABLE_CUSPARSE ON)
set(HYPRE_ENABLE_DEVICE_POOL OFF)
set(HYPRE_ENABLE_CUBLAS OFF)
set(HYPRE_ENABLE_CURAND ON)
set(HYPRE_ENABLE_GPU_PROFILING OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

if(NOT HYPRE_ENABLE_HYPRE_BLAS)
  find_dependency(BLAS)
endif()

if(NOT HYPRE_ENABLE_HYPRE_LAPACK)
  find_dependency(LAPACK)
endif()

if(HYPRE_USING_DSUPERLU)
  list(APPEND TPL_LIBRARIES OFF stdc++)
endif()

if(HYPRE_WITH_MPI)
  enable_language(C)
  find_dependency(MPI 3.1 EXACT COMPONENTS C)
endif()

if(HYPRE_WITH_OPENMP)
  find_dependency(OpenMP)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/HYPRETargets.cmake")
