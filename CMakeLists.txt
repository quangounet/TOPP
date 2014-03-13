# -*- coding: utf-8 -*-
# Copyright (C) 2014 Rosen Diankov <rosen.diankov@gmail.com> & Cuong Pham
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
cmake_minimum_required (VERSION 2.8.0)
project(topp)
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

set (TOPP_VERSION_MAJOR 0)
set (TOPP_VERSION_MINOR 1)
set (TOPP_VERSION_PATCH 0)
set (TOPP_VERSION ${TOPP_VERSION_MAJOR}.${TOPP_VERSION_MINOR}.${TOPP_VERSION_PATCH})
set (TOPP_SOVERSION ${TOPP_VERSION_MAJOR}.${TOPP_VERSION_MINOR})
set (CLIENT_SOVERSION ${TOPP_VERSION_MAJOR}.${TOPP_VERSION_MINOR})
message(STATUS "Compiling TOPP Version ${TOPP_VERSION}, soversion=${CLIENT_SOVERSION}")

set(TOPP_CMAKE_INSTALL_DIR "topp-${TOPP_VERSION_MAJOR}.${TOPP_VERSION_MINOR}" CACHE STRING "Directory to install the cmake config files.")
set(TOPP_TARGET_PROCESSOR ${CMAKE_SYSTEM_PROCESSOR} CACHE STRING "The target processor architecture to build for, this is combined with the generator toolchain")

message(STATUS "Using cmake version ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" )
# http://www.cmake.org/cmake/help/cmake-2.6.html#policy:CMP0002
cmake_policy(SET CMP0002 NEW)
# http://www.cmake.org/cmake/help/cmake-2.6.html#policy:CMP0003
cmake_policy(SET CMP0003 NEW)

# Use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# When building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(LIB_SUFFIX CACHE STRING "suffix for the library directory need for x86-64 systems that use lib64 ")

# The RPATH to be used when installing
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")

# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)
include(CheckFunctionExists)
include(CheckCXXSourceCompiles)
include(CheckCXXSourceRuns)
include(CheckCXXCompilerFlag)

# have to include before boost since the boost headers can be located in a previous installed version of this library
set(TOPP_INCLUDE_LOCAL_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(${TOPP_INCLUDE_LOCAL_DIRS})

if( NOT $ENV{BOOST_INCLUDEDIR} STREQUAL "" )
  set(Boost_INCLUDE_DIR $ENV{BOOST_INCLUDEDIR})
endif()
if( NOT $ENV{BOOST_LIBRARYDIR} STREQUAL "" )
  set(Boost_LIBRARY_DIRS $ENV{BOOST_LIBRARYDIR})
endif()
set(Boost_ADDITIONAL_VERSIONS "1.49" "1.48" "1.47" "1.46.1" "1.45" "1.44" "1.43" "1.42" "1.41" "1.40" "1.39" "1.38" "1.37.0" "1.37" "1.35.0" "1.34.1" "1.34.0" "1.34" "1.33.1" "1.33.0" "1.33" ${Boost_ADDITIONAL_VERSIONS})
#find_package(Boost COMPONENTS filesystem system thread date_time)

set(EXTRA_LIBRARIES ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
set(TOPP_LINK_DIRS)

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR COMPILER_IS_CLANG )
  set(CMAKE_CXX_FLAGS_OPTIMIZED "-O3 -DNDEBUG -DBOOST_DISABLE_ASSERTS -D_SECURE_SCL=0") # this practically removes all checks making it a very dangerous options to play with
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -D_DEBUG")
  
  add_definitions(-fPIC) # this is a very important switch and some libraries seem now to have it....
  # for some reason c++-0x is required...
  set(EXTRA_COMPILE_FLAGS "-fPIC -std=c++0x")
else()
  set(EXTRA_COMPILE_FLAGS "")
endif()

if( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
  add_definitions("-fno-strict-aliasing -Wall")
endif( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )

#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/topp/config.h IMMEDIATE @ONLY)

find_package(OpenRAVE)
if( OpenRAVE_FOUND )
  set(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -DWITH_OPENRAVE")
endif()

find_package(Boost ${OpenRAVE_Boost_VERSION} EXACT COMPONENTS iostreams python thread filesystem system thread date_time)# regex)
if( Boost_FOUND )
  include_directories(${Boost_INCLUDE_DIRS})
  set(TOPP_LINK_DIRS ${TOPP_LINK_DIRS} ${Boost_LIBRARY_DIRS})
elseif(Boost_VERSION AND NOT "${Boost_VERSION}" STREQUAL "0")
  include_directories(${Boost_INCLUDE_DIRS})
  set(TOPP_LINK_DIRS ${TOPP_LINK_DIRS} ${Boost_LIBRARY_DIRS})
  message(FATAL_ERROR "Could not find boost libraries!")
endif()

message(STATUS "found boost version: ${Boost_VERSION} at ${Boost_INCLUDE_DIRS}")
if( NOT Boost_THREAD_LIBRARY OR NOT Boost_DATE_TIME_LIBRARY )
  message(FATAL_ERROR "boost thread and date_time libraries are necessary")
endif()


include(FindPythonInterp)
## check python
find_package(PythonLibs) # using PYTHON_INCLUDE_PATH instead of PYTHON_INCLUDE_DIRS?
if( NOT PYTHON_EXECUTABLE )
  # look specifically for 2.6
  find_program(PYTHON_EXECUTABLE NAMES python2.6 python PATHS [HKEY_LOCAL_MACHINE\\SOFTWARE\\Python\\PythonCore\\2.6\\InstallPath])
endif()

# back compat?
if( NOT PYTHON_INCLUDE_DIRS )
  set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_PATH})
endif()

if( PYTHON_EXECUTABLE )
  get_filename_component(PYTHON_EXECUTABLE_NAME ${PYTHON_EXECUTABLE} NAME)
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_include()"
    OUTPUT_VARIABLE _python_numpy_include OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _python_failed0)
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; from os.path import relpath; print relpath(get_python_lib(1,prefix='${CMAKE_INSTALL_PREFIX}'),'${CMAKE_INSTALL_PREFIX}')"
    OUTPUT_VARIABLE TOPP_PYTHON_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    RESULT_VARIABLE _python_failed2)

  if( ${_python_failed0} EQUAL 0 AND ${_python_failed2} EQUAL 0 )
    if( NOT IS_DIRECTORY "${_python_numpy_include}" )
      set(PYTHON_EXECUTABLE)
    else()
      set(PYTHON_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS} ${_python_numpy_include})
    endif()

    # get the major.minor python version
    execute_process(
      COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print('%d.%d'%sys.version_info[0:2])"
      OUTPUT_VARIABLE _python_version OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _python_failed)
    if( ${_python_failed} EQUAL 0 )
      string(REGEX REPLACE "[\r\n]" "" PYTHON_MAJORMINOR_VERSION "${_python_version}")
    else()
      message(STATUS "failed to get python version")
    endif()
  else()
    message(STATUS "failed to get python site-package directories via get_python_lib")
    set(PYTHON_EXECUTABLE)
  endif()
endif()

set(TOPPPY_INSTALL_DIR "${TOPP_PYTHON_INSTALL_DIR}/TOPP" CACHE PATH "TOPP Python bindings installation directory")

find_package(PkgConfig) # pkg_check_modules

check_library_exists(lapack cgesdd_ "" LAPACK_FOUND)
if( LAPACK_FOUND )
  set(LAPACK_LIBRARIES lapack)
else()
  message(STATUS "lapack dev library not found")
endif()

add_subdirectory(src)

# add make uninstall capability
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" IMMEDIATE @ONLY)
#add_custom_target(uninstall "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/topp-config.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/topp-config.cmake" IMMEDIATE @ONLY)
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/topp-config-version.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/topp-config-version.cmake" IMMEDIATE @ONLY)

#install(FILES "${CMAKE_CURRENT_BINARY_DIR}/topp-config.cmake" "${CMAKE_CURRENT_BINARY_DIR}/topp-config-version.cmake" DESTINATION "lib${LIB_SUFFIX}/cmake/${TOPP_CMAKE_INSTALL_DIR}")
