/*  Copyright 2014 IST Austria

Contributed by: Jan Reininghaus

This file is part of DIPHA.

DIPHA is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DIPHA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with DIPHA.  If not, see <http://www.gnu.org/licenses/>. */


#define DIPHA_TEST
#include <string>
#include <gtest/gtest.h>
#include <random>
#include <limits>


#include <dipha/includes.h>




int main( int argc, char **argv )
{
    // initialize MPI environment
    MPI_Init( &argc, &argv );

    // redirect stdout to logfile - this is neccesary since we run multiple processes
    std::string filename = std::string( "unit_tests_rank_" ) + std::to_string( dipha::mpi_utils::get_rank( ) )
        + std::string( "_of_" ) + std::to_string( dipha::mpi_utils::get_num_processes( ) - 1 ) + std::string( ".log" );
    freopen( filename.c_str( ), "w", stdout );

    // initialize test framework
    testing::InitGoogleTest( &argc, argv );

    // run all tests
    int result = RUN_ALL_TESTS();

    // close logfile
    fclose( stdout );

    // shutdown MPI environment
    MPI_Finalize();

    // return test results
    return result;
}


