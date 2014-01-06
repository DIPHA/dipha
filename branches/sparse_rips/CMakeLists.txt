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

include_directories (include ${MPI_INCLUDE_PATH})

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

add_executable (test_distributed_sort tests/test_distributed_sort.cpp)
target_link_libraries(test_distributed_sort ${MPI_LIBRARIES})

source_group(Header\ Files FILES ${dipha_general_includes})
source_group(Header\ Files\\inputs FILES ${dipha_input_includes})
source_group(Header\ Files\\outputs FILES ${dipha_output_includes})
source_group(Header\ Files\\data_structures FILES ${dipha_data_structures_includes})
source_group(Header\ Files\\algorithms FILES ${dipha_algorithms_includes})

enable_testing()

set(NUM_PROCESSES 10)

add_test(NAME distributed_sort
         COMMAND mpiexec -n ${NUM_PROCESSES} $<TARGET_FILE:test_distributed_sort>
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )

add_test(NAME create_dualize
         COMMAND mpiexec -n ${NUM_PROCESSES} $<TARGET_FILE:dualize> primal_explicit.complex dualize.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )
add_test(NAME compare_dualize
         COMMAND ${CMAKE_COMMAND} -E compare_files dualize.DIPHA_TEST dual_explicit.complex
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )

add_test(NAME create_primal_cubical
         COMMAND mpiexec -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> cubical.complex primal_cubical.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )
add_test(NAME compare_primal_cubical
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_cubical.DIPHA_TEST cubical.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )

add_test(NAME create_dual_cubical
         COMMAND mpiexec -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual cubical.complex dual_cubical.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )
add_test(NAME compare_dual_cubical
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_cubical.DIPHA_TEST cubical.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )

add_test(NAME create_primal_explicit
         COMMAND mpiexec -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> primal_explicit.complex primal_explicit.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )
add_test(NAME compare_primal_explicit
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_explicit.DIPHA_TEST explicit.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )

add_test(NAME create_dual_explicit
         COMMAND mpiexec -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual dual_explicit.complex dual_explicit.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )
add_test(NAME compare_dual_explicit
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_explicit.DIPHA_TEST explicit.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )

add_test(NAME create_primal_extrinsic_full_rips
         COMMAND mpiexec -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> extrinsic_full_rips.complex primal_extrinsic_full_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )
add_test(NAME compare_primal_extrinsic_full_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_extrinsic_full_rips.DIPHA_TEST extrinsic_full_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )

add_test(NAME create_dual_extrinsic_full_rips
         COMMAND mpiexec -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual extrinsic_full_rips.complex dual_extrinsic_full_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )
add_test(NAME compare_dual_extrinsic_full_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_extrinsic_full_rips.DIPHA_TEST extrinsic_full_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )

add_test(NAME create_primal_intrinsic_full_rips
         COMMAND mpiexec -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> intrinsic_full_rips.complex primal_intrinsic_full_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )
add_test(NAME compare_primal_intrinsic_full_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_intrinsic_full_rips.DIPHA_TEST intrinsic_full_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )

add_test(NAME create_dual_intrinsic_full_rips
         COMMAND mpiexec -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual intrinsic_full_rips.complex dual_intrinsic_full_rips.DIPHA_TEST
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )
add_test(NAME compare_dual_intrinsic_full_rips
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_intrinsic_full_rips.DIPHA_TEST intrinsic_full_rips.diagram
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )