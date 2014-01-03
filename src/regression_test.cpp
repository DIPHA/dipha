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

    // Debugger attach point hack:
    //MPI_Barrier( MPI_COMM_WORLD );
    //if( dipha::mpi_utils::is_root() )
    //    system("pause");
    //MPI_Barrier( MPI_COMM_WORLD );

    // test distributed sorting
    {
        srand( (int)time( NULL ) + dipha::mpi_utils::get_rank() );
        long test_size = ( 1 << 13 ) - 1;
        std::vector< int > random_numbers( test_size );
        for( auto& random_number : random_numbers )
            random_number = rand();
        std::vector< long > distribution( dipha::mpi_utils::get_num_processes(), test_size );
        p_sort::parallel_sort( random_numbers.begin(), random_numbers.end(), distribution.data(), MPI_COMM_WORLD );
        std::vector< int > sorted_numbers( dipha::mpi_utils::get_num_processes() * test_size );
        MPI_Gather( random_numbers.data(), test_size, MPI_INT, sorted_numbers.data(), test_size, MPI_INT, 0, MPI_COMM_WORLD );
        if( dipha::mpi_utils::is_root() && !std::is_sorted( sorted_numbers.begin(), sorted_numbers.end() ) ) {
            dipha::mpi_utils::error_printer_if_root() << "Distributed sorting is buggy!" << std::endl;
            MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
        }
    }

    // test dualization of explicit complexes
    {
        dipha::algorithms::dualize( "primal_explicit.complex", "DIPHA_TEST_dual_explicit.complex" );
        MPI_Barrier( MPI_COMM_WORLD );
        bool result_has_not_changed = dipha::mpi_utils::are_files_equal( "DIPHA_TEST_dual_explicit.complex", "dual_explicit.complex" );
        MPI_Barrier( MPI_COMM_WORLD );
        
        if( dipha::mpi_utils::is_root() ) {
            std::remove( "DIPHA_TEST_dual_explicit.complex" );
            if( !result_has_not_changed ) {
                dipha::mpi_utils::error_printer_if_root() << "Dualization is buggy!" << std::endl;
                MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
            }
        }
    }

    // test primal reduction of cubical complex
    {
        dipha::inputs::weighted_cubical_complex complex;
        complex.load_binary( "cubical.complex" );
        dipha::data_structures::distributed_vector< int64_t > filtration_to_cell_map;
        dipha::data_structures::write_once_column_array reduced_columns;
        dipha::algorithms::compute_reduced_columns( complex, false, filtration_to_cell_map, reduced_columns );
        dipha::outputs::save_persistence_diagram( "DIPHA_TEST_cubical.diagram", complex, filtration_to_cell_map, reduced_columns, false );
        MPI_Barrier( MPI_COMM_WORLD );
        bool result_has_not_changed = dipha::mpi_utils::are_files_equal( "DIPHA_TEST_cubical.diagram", "cubical.diagram" );
        MPI_Barrier( MPI_COMM_WORLD );
        
        if( dipha::mpi_utils::is_root() ) {
            std::remove( "DIPHA_TEST_cubical.diagram" );
            if( !result_has_not_changed ) {
                dipha::mpi_utils::error_printer_if_root() << "Primal reduction of cubical_complex is buggy!" << std::endl;
                MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
            }
        }
    }

    // test dual reduction of cubical complex
    {
        dipha::inputs::weighted_cubical_complex complex;
        complex.load_binary( "cubical.complex" );
        dipha::data_structures::distributed_vector< int64_t > filtration_to_cell_map;
        dipha::data_structures::write_once_column_array reduced_columns;
        dipha::algorithms::compute_reduced_columns( complex, true, filtration_to_cell_map, reduced_columns );
        dipha::outputs::save_persistence_diagram( "DIPHA_TEST_cubical.diagram", complex, filtration_to_cell_map, reduced_columns, true );
        MPI_Barrier( MPI_COMM_WORLD );
        bool result_has_not_changed = dipha::mpi_utils::are_files_equal( "DIPHA_TEST_cubical.diagram", "cubical.diagram" );
        MPI_Barrier( MPI_COMM_WORLD );

        if( dipha::mpi_utils::is_root( ) ) {
            std::remove( "DIPHA_TEST_cubical.diagram" );
            if( !result_has_not_changed ) {
                dipha::mpi_utils::error_printer_if_root( ) << "Dual reduction of cubical_complex is buggy!" << std::endl;
                MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
            }
        }
    }

    // test primal reduction of explicit complex
    {
        dipha::inputs::weighted_explicit_complex complex; 
        complex.load_binary( "primal_explicit.complex" );
        dipha::data_structures::distributed_vector< int64_t > filtration_to_cell_map;
        dipha::data_structures::write_once_column_array reduced_columns;
        dipha::algorithms::compute_reduced_columns( complex, false, filtration_to_cell_map, reduced_columns );
        dipha::outputs::save_persistence_diagram( "DIPHA_TEST_explicit.diagram", complex, filtration_to_cell_map, reduced_columns, false );
        MPI_Barrier( MPI_COMM_WORLD );
        bool result_has_not_changed = dipha::mpi_utils::are_files_equal( "DIPHA_TEST_explicit.diagram", "explicit.diagram" );
        MPI_Barrier( MPI_COMM_WORLD );

        if( dipha::mpi_utils::is_root( ) ) {
            std::remove( "DIPHA_TEST_explicit.diagram" );
            if( !result_has_not_changed ) {
                dipha::mpi_utils::error_printer_if_root( ) << "Primal reduction of explicit_complex is buggy!" << std::endl;
                MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
            }
        }
    }

    // test dual reduction of explicit complex
    {
        dipha::inputs::weighted_explicit_complex complex;
        complex.load_binary( "dual_explicit.complex" );
        dipha::data_structures::distributed_vector< int64_t > filtration_to_cell_map;
        dipha::data_structures::write_once_column_array reduced_columns;
        dipha::algorithms::compute_reduced_columns( complex, true, filtration_to_cell_map, reduced_columns );
        dipha::outputs::save_persistence_diagram( "DIPHA_TEST_explicit.diagram", complex, filtration_to_cell_map, reduced_columns, true );
        MPI_Barrier( MPI_COMM_WORLD );
        bool result_has_not_changed = dipha::mpi_utils::are_files_equal( "DIPHA_TEST_explicit.diagram", "explicit.diagram" );
        MPI_Barrier( MPI_COMM_WORLD );

        if( dipha::mpi_utils::is_root( ) ) {
            std::remove( "DIPHA_TEST_explicit.diagram" );
            if( !result_has_not_changed ) {
                dipha::mpi_utils::error_printer_if_root( ) << "Dual reduction of explicit_complex is buggy!" << std::endl;
                MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
            }
        }
    }

    // test primal reduction of extrinsic full rips complex
    {
        dipha::inputs::full_rips_complex complex;
        complex.load_binary( "extrinsic_full_rips.complex" );
        dipha::data_structures::distributed_vector< int64_t > filtration_to_cell_map;
        dipha::data_structures::write_once_column_array reduced_columns;
        dipha::algorithms::compute_reduced_columns( complex, true, filtration_to_cell_map, reduced_columns );
        dipha::outputs::save_persistence_diagram( "DIPHA_TEST_extrinsic_full_rips.diagram", complex, filtration_to_cell_map, reduced_columns, true );
        MPI_Barrier( MPI_COMM_WORLD );
        bool result_has_not_changed = dipha::mpi_utils::are_files_equal( "DIPHA_TEST_extrinsic_full_rips.diagram", "extrinsic_full_rips.diagram" );
        MPI_Barrier( MPI_COMM_WORLD );

        if( dipha::mpi_utils::is_root( ) ) {
            std::remove( "DIPHA_TEST_extrinsic_full_rips.diagram" );
            if( !result_has_not_changed ) {
                dipha::mpi_utils::error_printer_if_root( ) << "Primal reduction of extrinsic_full_rips is buggy!" << std::endl;
                MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
            }
        }
    }

    // test dual reduction of extrinsic full rips complex
    {
        dipha::inputs::full_rips_complex complex;
        complex.load_binary( "extrinsic_full_rips.complex" );
        dipha::data_structures::distributed_vector< int64_t > filtration_to_cell_map;
        dipha::data_structures::write_once_column_array reduced_columns;
        dipha::algorithms::compute_reduced_columns( complex, false, filtration_to_cell_map, reduced_columns );
        dipha::outputs::save_persistence_diagram( "DIPHA_TEST_extrinsic_full_rips.diagram", complex, filtration_to_cell_map, reduced_columns, false );
        MPI_Barrier( MPI_COMM_WORLD );
        bool result_has_not_changed = dipha::mpi_utils::are_files_equal( "DIPHA_TEST_extrinsic_full_rips.diagram", "extrinsic_full_rips.diagram" );
        MPI_Barrier( MPI_COMM_WORLD );

        if( dipha::mpi_utils::is_root( ) ) {
            std::remove( "DIPHA_TEST_extrinsic_full_rips.diagram" );
            if( !result_has_not_changed ) {
                dipha::mpi_utils::error_printer_if_root( ) << "Dual reduction of extrinsic_full_rips is buggy!" << std::endl;
                MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
            }
        }
    }


    dipha::mpi_utils::cout_if_root() << std::endl << "All tests ran successfully." << std::endl;

    MPI_Finalize();
}
