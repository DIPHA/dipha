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

void print_help()
{
    std::cerr << "Usage: " << "dipha " << "[options] input_filename output_filename" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << std::endl;
    std::cerr << "--help    --  prints this screen" << std::endl;
    std::cerr << "--dual    --  use dualization" << std::endl;
    std::cerr << "--benchmark --  prints timing info" << std::endl;
}

void print_help_and_exit()
{
    print_help();
    MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
}

void parse_command_line( int argc, char** argv, bool& benchmark, bool& dualize, std::string& input_filename, std::string& output_filename )
{

    if( argc < 3 )
        print_help_and_exit();

    input_filename = argv[ argc - 2 ];
    output_filename = argv[ argc - 1 ];

    for( int idx = 1; idx < argc - 2; idx++ ) {
        const std::string option = argv[ idx ];
        if( option == "--benchmark" ) benchmark = true;
        else if( option == "--help" ) print_help_and_exit();
        else if( option == "--dual" ) dualize = true;
        else print_help_and_exit();
    }
}

template< typename Complex >
void compute( const std::string& input_filename,
              bool dualize,
              const std::string& output_filename )
{
    Complex complex;
    DIPHA_MACROS_BENCHMARK( complex.load_binary( input_filename ); );
    dipha::data_structures::distributed_vector< int64_t > filtration_to_cell_map;
    dipha::data_structures::write_once_column_array reduced_columns;
    dipha::algorithms::compute_reduced_columns( complex, dualize, filtration_to_cell_map, reduced_columns );
    DIPHA_MACROS_BENCHMARK( dipha::outputs::save_persistence_diagram( output_filename, complex, filtration_to_cell_map, reduced_columns, dualize ); );
}

int main( int argc, char** argv )
{
    // mandatory MPI initilization call
    MPI_Init( &argc, &argv );

    // take time at beggining of execution
    double time_at_start = MPI_Wtime();

    std::string input_filename; // name of file that contains the weighted cell complex
    std::string output_filename; // name of file that will contain the persistence diagram
    bool benchmark = false; // print timings / info
    bool dualize = false; // primal / dual computation toggle
    parse_command_line( argc, argv, benchmark, dualize, input_filename, output_filename );

    if( benchmark ) {
        dipha::globals::benchmark = true;

        dipha::mpi_utils::cout_if_root() << std::endl << "Input filename: " << std::endl << input_filename << std::endl;
        dipha::mpi_utils::cout_if_root() << std::endl << "Number of processes used: " << std::endl << dipha::mpi_utils::get_num_processes() << std::endl;
        dipha::mpi_utils::cout_if_root() << std::endl << "Detailed information for rank 0:" << std::endl;
        dipha::mpi_utils::cout_if_root() << std::setw( 11 ) << "time" << std::setw( 13 ) << "prior mem" << std::setw( 13 ) << "peak mem" << std::setw( 13 ) << "bytes recv" << std::endl;
    }

    switch( dipha::file_types::get_file_type( input_filename ) ) {
    case dipha::file_types::WEIGHTED_CUBICAL_COMPLEX:
        compute< dipha::inputs::weighted_cubical_complex >( input_filename, dualize, output_filename );
        break;
    case dipha::file_types::WEIGHTED_EXPLICIT_COMPLEX:
        compute< dipha::inputs::weighted_explicit_complex >( input_filename, dualize, output_filename );
        break;
    case dipha::file_types::EXTRINSIC_FULL_RIPS_COMPLEX:
      // Go to next case
    case dipha::file_types::INTRINSIC_FULL_RIPS_COMPLEX:
        compute< dipha::inputs::full_rips_complex >( input_filename, dualize, output_filename );
        break;
    default:
        dipha::mpi_utils::error_printer_if_root() << "Unknown complex type in DIPHA file" << input_filename << std::endl;
        MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
    }

    if( benchmark ) {
        MPI_Barrier( MPI_COMM_WORLD );
        dipha::mpi_utils::cout_if_root() << std::endl << "Overall running time in seconds: " << std::endl;
        dipha::mpi_utils::cout_if_root() << std::setprecision( 1 ) << MPI_Wtime() - time_at_start << std::endl;

        dipha::mpi_utils::cout_if_root() << std::endl << "Reduction kernel running time in seconds: " << std::endl;
        dipha::mpi_utils::cout_if_root() << std::setprecision( 1 ) << dipha::globals::reduction_kernel_running_time << std::endl;

        int64_t peak_mem = getPeakRSS() >> 20;
        std::vector< int64_t > peak_mem_per_rank( dipha::mpi_utils::get_num_processes() );
        MPI_Gather( &peak_mem, 1, MPI_LONG_LONG, peak_mem_per_rank.data(), 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
        dipha::mpi_utils::cout_if_root() << std::endl << "Overall peak mem in GB of all ranks: " << std::endl;
        dipha::mpi_utils::cout_if_root() << (double)*std::max_element( peak_mem_per_rank.begin(), peak_mem_per_rank.end() ) / 1024.0 << std::endl;

        std::vector< int64_t > bytes_received_per_rank( dipha::mpi_utils::get_num_processes() );
        MPI_Gather( &dipha::globals::bytes_received, 1, MPI_LONG_LONG, bytes_received_per_rank.data(), 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );

        dipha::mpi_utils::cout_if_root() << std::endl << "Maximal communication traffic (without sorting) in GB between any pair of nodes:" << std::endl;
        dipha::mpi_utils::cout_if_root() << std::setprecision( 1 ) << (double)( *std::max_element( bytes_received_per_rank.begin(), bytes_received_per_rank.end() ) >> 20 ) / 1024.0 << std::endl;

        dipha::mpi_utils::cout_if_root() << std::endl << "Total communication traffic (without sorting) in GB between all pairs of nodes:" << std::endl;
        dipha::mpi_utils::cout_if_root() << std::setprecision( 1 ) << (double)( std::accumulate( bytes_received_per_rank.begin(), bytes_received_per_rank.end(), 0LL ) >> 20 ) / 1024.0 << std::endl;
    }

    MPI_Finalize();
}
