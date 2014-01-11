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

TEST( distributed_sort, psort )
{
    std::minstd_rand generator( ::testing::UnitTest::GetInstance( )->random_seed( ) );
    std::uniform_int_distribution< int > distribution( 0, 1 << 16 );
    long test_size = distribution( generator );

    std::vector< int > random_numbers( test_size );
    for( auto& random_number : random_numbers )
        random_number = distribution( generator );
    std::vector< long > elem_distribution( mpi_utils::get_num_processes( ), test_size );
    p_sort::parallel_sort( random_numbers.begin( ), random_numbers.end( ), elem_distribution.data( ), MPI_COMM_WORLD );
    std::vector< int > sorted_numbers( mpi_utils::get_num_processes( ) * test_size );
    MPI_Gather( random_numbers.data( ), test_size, MPI_INT, sorted_numbers.data( ), test_size, MPI_INT, 0, MPI_COMM_WORLD );

    if( mpi_utils::is_root( ) )
        ASSERT_TRUE( std::is_sorted( sorted_numbers.begin( ), sorted_numbers.end( ) ) );
}