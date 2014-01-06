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

#include <dipha/includes.h> 

int main( int argc, char** argv )
{
    MPI_Init( &argc, &argv );
    
    srand( (int)time( NULL ) + dipha::mpi_utils::get_rank() );
    long test_size = ( 1 << 13 ) - 1;
    std::vector< int > random_numbers( test_size );
    for( auto& random_number : random_numbers )
        random_number = rand();
    std::vector< long > distribution( dipha::mpi_utils::get_num_processes(), test_size );
    p_sort::parallel_sort( random_numbers.begin(), random_numbers.end(), distribution.data(), MPI_COMM_WORLD );
    std::vector< int > sorted_numbers( dipha::mpi_utils::get_num_processes() * test_size );
    MPI_Gather( random_numbers.data(), test_size, MPI_INT, sorted_numbers.data(), test_size, MPI_INT, 0, MPI_COMM_WORLD );
    
    if( dipha::mpi_utils::is_root() && !std::is_sorted( sorted_numbers.begin(), sorted_numbers.end() ) )
        MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );

    MPI_Finalize();
}
