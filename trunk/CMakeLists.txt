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
    
    # we need the following for gtest
    if( MSVC_VERSION EQUAL 1700 )
        add_definitions(-D_VARIADIC_MAX=10)
    endif() 

	# gtest need static runtime lib
    foreach(flag_var
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
        if(${flag_var} MATCHES "/MD")
            string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
        endif(${flag_var} MATCHES "/MD")
    endforeach(flag_var)
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    add_definitions(-std=c++11 -w)
endif()

FIND_PACKAGE(MPI REQUIRED)

add_subdirectory(externals/gtest-1.7.0)

include_directories (include externals ${MPI_INCLUDE_PATH} ${gtest_SOURCE_DIR}/include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${GTEST_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS} ${GTEST_LINK_FLAGS}")

FILE(GLOB_RECURSE all_dipha_includes "include/dipha/*.h")
FILE(GLOB dipha_general_includes "include/dipha/*.h")
FILE(GLOB dipha_data_structures_includes "include/dipha/data_structures/*.h")
FILE(GLOB dipha_input_includes "include/dipha/inputs/*.h")
FILE(GLOB dipha_output_includes "include/dipha/outputs/*.h")
FILE(GLOB dipha_algorithms_includes "include/dipha/algorithms/*.h")

add_executable (dipha src/dipha.cpp ${all_dipha_includes})
target_link_libraries(dipha gtest ${MPI_LIBRARIES})

add_executable (dualize src/dualize.cpp)
target_link_libraries(dualize gtest ${MPI_LIBRARIES})

add_executable(unit_tests tests/unit_tests.cpp)
target_link_libraries(unit_tests gtest ${MPI_LIBRARIES})

source_group(Header\ Files FILES ${dipha_general_includes})
source_group(Header\ Files\\inputs FILES ${dipha_input_includes})
source_group(Header\ Files\\outputs FILES ${dipha_output_includes})
source_group(Header\ Files\\data_structures FILES ${dipha_data_structures_includes})
source_group(Header\ Files\\algorithms FILES ${dipha_algorithms_includes})

enable_testing()

string(RANDOM LENGTH 4 ALPHABET 0123456789 RANDOM_SEED)
math(EXPR RANDOM_SEED "${RANDOM_SEED} + 1")  
message( STATUS "RANDOM_SEED = ${RANDOM_SEED}") 

string(RANDOM LENGTH 1 ALPHABET 0123456789 NUM_PROCESSES)
math(EXPR NUM_PROCESSES "${NUM_PROCESSES} + 1") 
message( STATUS "NUM_PROCESSES = ${NUM_PROCESSES}") 

set(NUM_TEST_RUNS 10)
message( STATUS "NUM_TEST_RUNS = ${NUM_TEST_RUNS}") 

FILE(GLOB unit_test_log_files "test_data/*.log")
if(NOT( unit_test_log_files STREQUAL "" ) )
	execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${unit_test_log_files}) 
endif()

FILE(GLOB functional_test_temp_data "test_data/*.DIPHA_TEST")
if(NOT( functional_test_temp_data STREQUAL "" ) )
	execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${functional_test_temp_data}) 
endif()

add_test(NAME unit_tests
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:unit_tests> --gtest_shuffle 
						 --gtest_repeat=${NUM_TEST_RUNS} --gtest_break_on_failure --gtest_random_seed=${RANDOM_SEED}
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )

add_test(NAME create_dualize
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dualize> primal_explicit.complex dualize.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_dualize
         COMMAND ${CMAKE_COMMAND} -E compare_files dualize.DIPHA_TEST dual_explicit.complex
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )

add_test(NAME create_primal_cubical
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> cubical.complex primal_cubical.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_primal_cubical
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_cubical.DIPHA_TEST cubical.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )

add_test(NAME create_dual_cubical
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual cubical.complex dual_cubical.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_dual_cubical
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_cubical.DIPHA_TEST cubical.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )

add_test(NAME create_primal_explicit
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> primal_explicit.complex primal_explicit.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_primal_explicit
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_explicit.DIPHA_TEST explicit.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )

add_test(NAME create_dual_explicit
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual dual_explicit.complex dual_explicit.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_dual_explicit
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_explicit.DIPHA_TEST explicit.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )

add_test(NAME create_primal_extrinsic_full_rips
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> extrinsic_full_rips.complex primal_extrinsic_full_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_primal_extrinsic_full_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_extrinsic_full_rips.DIPHA_TEST extrinsic_full_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )

add_test(NAME create_dual_extrinsic_full_rips
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual extrinsic_full_rips.complex dual_extrinsic_full_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_dual_extrinsic_full_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_extrinsic_full_rips.DIPHA_TEST extrinsic_full_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )

add_test(NAME create_primal_intrinsic_full_rips
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> intrinsic_full_rips.complex primal_intrinsic_full_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_primal_intrinsic_full_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_intrinsic_full_rips.DIPHA_TEST intrinsic_full_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )

add_test(NAME create_dual_intrinsic_full_rips
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual intrinsic_full_rips.complex dual_intrinsic_full_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )
add_test(NAME compare_dual_intrinsic_full_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_intrinsic_full_rips.DIPHA_TEST intrinsic_full_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/test_data
        )