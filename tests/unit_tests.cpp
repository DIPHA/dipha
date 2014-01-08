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

#include <gtest/gtest.h>

#include <dipha/includes.h>

#include <random>
#include <limits>

using namespace dipha;

TEST( distributed_sort, psort )
{
    std::minstd_rand generator( ::testing::UnitTest::GetInstance( )->random_seed( ) );
    std::uniform_int_distribution< int > distribution( 0, 1 << 16 );
    long test_size = distribution( generator );

    std::vector< int > random_numbers( test_size );
    for( auto& random_number : random_numbers )
        random_number = distribution( generator );
    std::vector< long > elem_distribution( mpi_utils::get_num_processes(), test_size );
    p_sort::parallel_sort( random_numbers.begin( ), random_numbers.end( ), elem_distribution.data( ), MPI_COMM_WORLD );
    std::vector< int > sorted_numbers( mpi_utils::get_num_processes() * test_size );
    MPI_Gather( random_numbers.data(), test_size, MPI_INT, sorted_numbers.data(), test_size, MPI_INT, 0, MPI_COMM_WORLD );

    if( mpi_utils::is_root() )
        ASSERT_TRUE( std::is_sorted( sorted_numbers.begin(), sorted_numbers.end() ) );
}


TEST( mpi_utils, send_receive_vector )
{
    std::minstd_rand generator( ::testing::UnitTest::GetInstance( )->random_seed( ) );
    std::uniform_int_distribution< int > distribution( 0, 1 << 16 );
    long test_size = distribution( generator );
    
    std::vector< int > message_to_send( test_size, mpi_utils::get_rank( ) );
    std::vector< MPI_Request > requests;
    for( int target_rank = 0; target_rank < mpi_utils::get_num_processes( ); target_rank++ ) {
        mpi_utils::non_blocking_send_vector( message_to_send, target_rank, mpi_utils::MSG_TESTING, requests );
    }

    std::vector< int > receive_buffer;
    for( int target_rank = 0; target_rank < mpi_utils::get_num_processes(); target_rank++ ) {
        mpi_utils::receive_vector( receive_buffer, target_rank, mpi_utils::MSG_TESTING );
        ASSERT_EQ( test_size, receive_buffer.size( ) );
        if( !receive_buffer.empty() ) {
            ASSERT_EQ( target_rank, *std::min_element( receive_buffer.cbegin( ), receive_buffer.cend( ) ) );
            ASSERT_EQ( target_rank, *std::max_element( receive_buffer.cbegin( ), receive_buffer.cend( ) ) );
        }
    }

    MPI_Waitall( (int)requests.size(), requests.data(), MPI_STATUSES_IGNORE );
}



int main( int argc, char **argv )
{
    // initialize MPI environment
    MPI_Init( &argc, &argv );

    // redirect stdout to logfile - this is neccesary since we run multiple processes
    std::string filename = std::string( "unit_tests_rank_" ) + std::to_string( mpi_utils::get_rank( ) )
                         + std::string( "_of_" ) + std::to_string( mpi_utils::get_num_processes( ) - 1 ) + std::string( ".log" );
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


