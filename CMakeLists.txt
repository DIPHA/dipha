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

add_executable (create_phat_filtration src/create_phat_filtration.cpp)
target_link_libraries(create_phat_filtration ${MPI_LIBRARIES})

add_executable (full_to_sparse_distance_matrix src/full_to_sparse_distance_matrix.cpp)
target_link_libraries(full_to_sparse_distance_matrix ${MPI_LIBRARIES})

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

SET(TEST_DATA_DIR ${CMAKE_HOME_DIRECTORY}/test_data)

macro(REGRESSION_TEST TEST_NAME UPPER_DIM INPUT)
    add_test(NAME create_primal_${TEST_NAME}
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --upper_dim ${UPPER_DIM} ${INPUT}.complex primal_${TEST_NAME}.DIPHA_TEST
         WORKING_DIRECTORY ${TEST_DATA_DIR}
        )
	add_test(NAME compare_primal_${TEST_NAME}
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_${TEST_NAME}.DIPHA_TEST ${TEST_NAME}.diagram
         WORKING_DIRECTORY ${TEST_DATA_DIR}
        )
	set_tests_properties(compare_primal_${TEST_NAME} PROPERTIES DEPENDS create_primal_${TEST_NAME})

	add_test(NAME create_dual_${TEST_NAME}
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual --upper_dim ${UPPER_DIM} ${INPUT}.complex dual_${TEST_NAME}.DIPHA_TEST
         WORKING_DIRECTORY ${TEST_DATA_DIR}
        )
	add_test(NAME compare_dual_${TEST_NAME}
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_${TEST_NAME}.DIPHA_TEST ${TEST_NAME}.diagram
         WORKING_DIRECTORY ${TEST_DATA_DIR}
        )
	set_tests_properties(compare_dual_${TEST_NAME} PROPERTIES DEPENDS create_dual_${TEST_NAME})
endmacro(REGRESSION_TEST)

SET(DEFAULT_UPPER_DIM 9223372036854775807)

REGRESSION_TEST(smooth_def ${DEFAULT_UPPER_DIM} smooth_16)     
REGRESSION_TEST(smooth_2 2 smooth_16)  

REGRESSION_TEST(noise_3_16_def ${DEFAULT_UPPER_DIM} noise_3_16)  
REGRESSION_TEST(noise_4_8_def ${DEFAULT_UPPER_DIM} noise_4_8)  

REGRESSION_TEST(ramp_3_16_def ${DEFAULT_UPPER_DIM} ramp_3_16)  
REGRESSION_TEST(ramp_4_8_def ${DEFAULT_UPPER_DIM} ramp_4_8)  

REGRESSION_TEST(sphere_3 3 sphere_3_32)
REGRESSION_TEST(sphere_2 2 sphere_3_32)
REGRESSION_TEST(sphere_1 1 sphere_3_32)

REGRESSION_TEST(sphere_3_2   3 sphere_3_32_2)  
REGRESSION_TEST(sphere_3_0.8 3 sphere_3_32_0.8)
REGRESSION_TEST(sphere_2_0.5 2 sphere_3_32_0.5)
REGRESSION_TEST(sphere_1_0.3 2 sphere_3_32_0.3)


add_test(NAME create_dualize
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dualize> primal_explicit.complex dualize.DIPHA_TEST
         WORKING_DIRECTORY ${TEST_DATA_DIR}
        )
add_test(NAME compare_dualize
         COMMAND ${CMAKE_COMMAND} -E compare_files dualize.DIPHA_TEST dual_explicit.complex
         WORKING_DIRECTORY ${TEST_DATA_DIR}
        )
set_tests_properties(compare_dualize PROPERTIES DEPENDS create_dualize)

add_test(NAME create_primal_explicit
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> primal_explicit.complex primal_explicit.DIPHA_TEST
         WORKING_DIRECTORY ${TEST_DATA_DIR}
        )
add_test(NAME compare_primal_explicit
         COMMAND ${CMAKE_COMMAND} -E compare_files primal_explicit.DIPHA_TEST explicit.diagram
         WORKING_DIRECTORY ${TEST_DATA_DIR}
        )
set_tests_properties(compare_primal_explicit PROPERTIES DEPENDS create_primal_explicit)

add_test(NAME create_dual_explicit
         COMMAND ${MPIEXEC} -n ${NUM_PROCESSES} $<TARGET_FILE:dipha> --dual dual_explicit.complex dual_explicit.DIPHA_TEST
         WORKING_DIRECTORY ${TEST_DATA_DIR}
        )
add_test(NAME compare_dual_explicit
         COMMAND ${CMAKE_COMMAND} -E compare_files dual_explicit.DIPHA_TEST explicit.diagram
         WORKING_DIRECTORY ${TEST_DATA_DIR}
        )
set_tests_properties(compare_dual_explicit PROPERTIES DEPENDS create_dual_explicit)      
