# Copyright (C) |Meso|Star> 2015-2018 (contact@meso-star.com)
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use,
# modify and/or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

cmake_minimum_required(VERSION 2.8)
project(model_problem)
enable_language(C)
set(MODEL_PROBLEM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/../src/)

################################################################################
# Check dependencies
################################################################################
find_package(RCMake 0.3 REQUIRED)
find_package(RSys 0.6 REQUIRED)
find_package(Star3D 0.4.1 REQUIRED)
find_package(Star3DAW 0.1.3 REQUIRED)
find_package(StarSP 0.7 REQUIRED)
find_package(StarMC 0.4 REQUIRED)

include_directories(
  ${RSys_INCLUDE_DIR}
  ${Star3D_INCLUDE_DIR}
  ${Star3DAW_INCLUDE_DIR}
  ${StarSP_INCLUDE_DIR}
  ${StarMC_INCLUDE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${RCMAKE_SOURCE_DIR})
include(rcmake)
include(rcmake_runtime)

################################################################################
# Configure and define the targets
################################################################################
set(VERSION_MAJOR 0)
set(VERSION_MINOR 3)
set(VERSION_PATCH 0)
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

set(MODEL_PROBLEM_FILES_SRC model_problem.c model_problem_realization.c)
set(MODEL_PROBLEM_FILES_INC model_problem_realization.h)
set(MODEL_PROBLEM_FILES_DOC COPYING.fr COPYING.en README.md)

# Prepend each file in the `MODEL_PROBLEM_FILES_<SRC|INC>' list by the absolute
# path of the directory in which they lie
rcmake_prepend_path(MODEL_PROBLEM_FILES_SRC ${MODEL_PROBLEM_SOURCE_DIR})
rcmake_prepend_path(MODEL_PROBLEM_FILES_INC ${MODEL_PROBLEM_SOURCE_DIR})
rcmake_prepend_path(MODEL_PROBLEM_FILES_DOC ${PROJECT_SOURCE_DIR}/../)

if(CMAKE_COMPILER_IS_GNUCC)
  set(MATH_LIB m)
endif()

add_executable(model_problem ${MODEL_PROBLEM_FILES_SRC} ${MODEL_PROBLEM_FILES_INC})
target_link_libraries(model_problem RSys Star3D Star3DAW StarSP StarMC ${MATH_LIB})
set_target_properties(model_problem PROPERTIES
  VERSION ${VERSION}
  SOVERSION ${VERSION_MAJOR})

################################################################################
# Define output & install directories
################################################################################
install(TARGETS model_problem
  ARCHIVE DESTINATION bin
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)
install(FILES ${MODEL_PROBLEM_FILES_DOC} DESTINATION share/doc/star-model_problem)

