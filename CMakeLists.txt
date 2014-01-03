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

add_executable (regression_test src/regression_test.cpp ${all_dipha_includes})
target_link_libraries(regression_test ${MPI_LIBRARIES})
add_executable (dualize src/dualize.cpp ${all_dipha_includes})
target_link_libraries(dualize ${MPI_LIBRARIES})
add_executable (dipha src/dipha.cpp ${all_dipha_includes})
target_link_libraries(dipha ${MPI_LIBRARIES})

source_group(Header\ Files FILES ${dipha_general_includes})
source_group(Header\ Files\\inputs FILES ${dipha_input_includes})
source_group(Header\ Files\\outputs FILES ${dipha_output_includes})
source_group(Header\ Files\\data_structures FILES ${dipha_data_structures_includes})
source_group(Header\ Files\\algorithms FILES ${dipha_algorithms_includes})

enable_testing()
get_target_property(regression_test_path regression_test LOCATION)
add_test(NAME regression_test
         COMMAND mpiexec -n 3 ${regression_test_path}
         WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}/regression_data
        )