#  Copyright 2014 IST Austria
#
#  Contributed by: Ulrich Bauer, Jan Reininghaus
#
#  This file is part of DIPHA.
#
#  DIPHA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  DIPHA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with DIPHA.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)

project("Dipha")

set(CMAKE_BUILD_TYPE "Release")

if(MSVC OR MSVC_IDE)
    if( MSVC_VERSION LESS 1700 )
        message(FATAL_ERROR "The project requires C++11 features. " 
        "You need at least Visual Studio 11 (Microsoft Visual Studio 2012)")
    endif()        
    
    # Compilation using Microsoft Visual Studio 2013 seems to be stuck in an infinite loop ...
    if( MSVC_VERSION GREATER 1700 )
        set(CMAKE_GENERATOR_TOOLSET "v110" CACHE STRING "Platform Toolset" FORCE)               
    endif()
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    add_definitions(-std=c++11 -w)
endif()

FIND_PACKAGE(MPI REQUIRED)

include_directories (include externals ${MPI_INCLUDE_PATH})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")

FILE(GLOB_RECURSE all_dipha_includes "include/dipha/*.h")
FILE(GLOB dipha_general_includes "include/dipha/*.h")
FILE(GLOB dipha_data_structures_includes "include/dipha/data_structures/*.h")
FILE(GLOB dipha_input_includes "include/dipha/inputs/*.h")
FILE(GLOB dipha_output_includes "include/dipha/outputs/*.h")
FILE(GLOB dipha_algorithms_includes "include/dipha/algorithms/*.h")

add_executable (dipha src/dipha.cpp ${all_dipha_includes})
target_link_libraries(dipha ${MPI_LIBRARIES})

add_executable (dualize src/dualize.cpp)
target_link_libraries(dualize ${MPI_LIBRARIES})

source_group(Header\ Files FILES ${dipha_general_includes})
source_group(Header\ Files\\inputs FILES ${dipha_input_includes})
source_group(Header\ Files\\outputs FILES ${dipha_output_includes})
source_group(Header\ Files\\data_structures FILES ${dipha_data_structures_includes})
source_group(Header\ Files\\algorithms FILES ${dipha_algorithms_includes})

enable_testing()

string(RANDOM LENGTH 1 ALPHABET 0123456789 NUM_PROCESSES)
math(EXPR NUM_PROCESSES "${NUM_PROCESSES} + 1") 
#set(NUM_PROCESSES 1)
message( STATUS "NUM_PROCESSES = ${NUM_PROCESSES}") 

FILE(GLOB functional_test_temp_data "test_data/*.DIPHA_TEST")
if(NOT( functional_test_temp_data STREQUAL "" ) )
	execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${functional_test_temp_data}) 
endif()

add_test(NAME create_dualize
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dualize> primal_explicit.complex dualize.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_dualize
         COMMAND ${CMAKE_COMMAND} -E compare_files dualize.DIPHA_TEST dual_explicit.complex
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
set_tests_properties(compare_dualize PROPERTIES DEPENDS create_dualize)

add_test(NAME create_primal_cubical
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> cubical.complex primal_cubical.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_primal_cubical
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_cubical.DIPHA_TEST cubical.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
set_tests_properties(compare_primal_cubical PROPERTIES DEPENDS create_primal_cubical)

add_test(NAME create_dual_cubical
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual cubical.complex dual_cubical.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_dual_cubical
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_cubical.DIPHA_TEST cubical.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
set_tests_properties(compare_dual_cubical PROPERTIES DEPENDS create_dual_cubical)
        
add_test(NAME create_primal_explicit
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> primal_explicit.complex primal_explicit.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_primal_explicit
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_explicit.DIPHA_TEST explicit.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
set_tests_properties(compare_primal_explicit PROPERTIES DEPENDS create_primal_explicit)

add_test(NAME create_dual_explicit
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual dual_explicit.complex dual_explicit.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_dual_explicit
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_explicit.DIPHA_TEST explicit.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
set_tests_properties(compare_dual_explicit PROPERTIES DEPENDS create_dual_explicit)

add_test(NAME create_primal_intrinsic_full_rips
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --upper_dim 2 intrinsic_full_rips.complex primal_intrinsic_full_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_primal_intrinsic_full_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_intrinsic_full_rips.DIPHA_TEST intrinsic_full_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
set_tests_properties(compare_primal_intrinsic_full_rips PROPERTIES DEPENDS create_primal_intrinsic_full_rips)           

add_test(NAME create_dual_intrinsic_full_rips
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual --upper_dim 2 intrinsic_full_rips.complex dual_intrinsic_full_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_dual_intrinsic_full_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_intrinsic_full_rips.DIPHA_TEST intrinsic_full_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
set_tests_properties(compare_dual_intrinsic_full_rips PROPERTIES DEPENDS create_dual_intrinsic_full_rips)           

add_test(NAME create_primal_intrinsic_sparse_rips_large_threshold
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> intrinsic_sparse_rips_large_threshold.complex primal_intrinsic_sparse_rips_large_threshold.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_primal_intrinsic_sparse_rips_large_threshold
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_intrinsic_sparse_rips_large_threshold.DIPHA_TEST intrinsic_sparse_rips_large_threshold.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
set_tests_properties(compare_primal_intrinsic_sparse_rips_large_threshold PROPERTIES DEPENDS create_primal_intrinsic_sparse_rips_large_threshold)             

add_test(NAME create_dual_intrinsic_sparse_rips_large_threshold
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual intrinsic_sparse_rips_large_threshold.complex dual_intrinsic_sparse_rips_large_threshold.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_dual_intrinsic_sparse_rips_large_threshold
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_intrinsic_sparse_rips_large_threshold.DIPHA_TEST intrinsic_sparse_rips_large_threshold.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
set_tests_properties(compare_dual_intrinsic_sparse_rips_large_threshold PROPERTIES DEPENDS create_dual_intrinsic_sparse_rips_large_threshold)              

add_test(NAME create_primal_intrinsic_sparse_rips
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> intrinsic_sparse_rips.complex primal_intrinsic_sparse_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_primal_intrinsic_sparse_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_intrinsic_sparse_rips.DIPHA_TEST intrinsic_sparse_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
set_tests_properties(compare_primal_intrinsic_sparse_rips PROPERTIES DEPENDS create_primal_intrinsic_sparse_rips)             

add_test(NAME create_dual_intrinsic_sparse_rips
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual intrinsic_sparse_rips.complex dual_intrinsic_sparse_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_dual_intrinsic_sparse_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_intrinsic_sparse_rips.DIPHA_TEST intrinsic_sparse_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
set_tests_properties(compare_dual_intrinsic_sparse_rips PROPERTIES DEPENDS create_dual_intrinsic_sparse_rips) 