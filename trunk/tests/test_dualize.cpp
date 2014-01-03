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

    dipha::algorithms::dualize( "primal_explicit.complex", "DIPHA_TEST_dual_explicit.complex" );
    MPI_Barrier( MPI_COMM_WORLD );
    bool result_has_not_changed = dipha::mpi_utils::are_files_equal( "DIPHA_TEST_dual_explicit.complex", "dual_explicit.complex" );
    MPI_Barrier( MPI_COMM_WORLD );
        
    if( dipha::mpi_utils::is_root() ) {
        std::remove( "DIPHA_TEST_dual_explicit.complex" );
        if( !result_has_not_changed ) {
            MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
        }
    }

    MPI_Finalize();
}
