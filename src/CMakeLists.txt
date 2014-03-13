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

# copy the python files

# include python bindings sample
if( Boost_PYTHON_FOUND AND Boost_THREAD_FOUND )
  if( PYTHONLIBS_FOUND OR PYTHON_LIBRARIES )
    if( PYTHON_EXECUTABLE )
      # get the site-packages directory
      execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print get_python_lib(1)"
        OUTPUT_VARIABLE _python_sitepackage
        RESULT_VARIABLE _python_failed)
      if( ${_python_failed} EQUAL 0 )
        string(REGEX REPLACE "[\r\n]" "" _python_sitepackage "${_python_sitepackage}")
        set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_PATH} ${_python_sitepackage}/numpy/core/include)
      else()
        message(STATUS "failed to get python site-package directory")
      endif()
    endif()
    
    include_directories(${PYTHON_INCLUDE_PATH} ${OpenRAVE_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
    link_directories(${OpenRAVE_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})
    file(GLOB topp_files ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
    
    add_library(TOPPbindings SHARED ${topp_files})
    target_link_libraries(TOPPbindings ${OpenRAVE_PYTHON_DIR}/openravepy_int.so ${PYTHON_LIBRARIES} ${Boost_PYTHON_LIBRARY} ${Boost_THREAD_LIBRARY} ${OpenRAVE_LIBRARIES})
    set_target_properties(TOPPbindings PROPERTIES PREFIX "" COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} ${OpenRAVE_CXX_FLAGS}" LINK_FLAGS "${OpenRAVE_LINK_FLAGS}")
    if( WIN32 )
      set_target_properties(TOPPbindings PROPERTIES SUFFIX ".pyd")
    endif()
    install(TARGETS TOPPbindings DESTINATION ${TOPPPY_INSTALL_DIR} COMPONENT ${COMPONENT_PREFIX}python)
  endif()
endif()

file(GLOB topp_python_files "${CMAKE_CURRENT_SOURCE_DIR}/python/*.py")
install(FILES ${topp_python_files} DESTINATION ${TOPPPY_INSTALL_DIR} COMPONENT ${COMPONENT_PREFIX}python)

#install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python/" DESTINATION ${TOPPPY_INSTALL_DIR} FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ COMPONENT ${COMPONENT_PREFIX}python PATTERN ".svn" EXCLUDE PATTERN ".pyc" EXCLUDE)
