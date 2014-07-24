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

#pragma once

#include <dipha/includes.h>

namespace dipha {
    namespace outputs {
        /** file format: file_types::DIPHA % file_types::PERSISTENCE_DIAGRAM % global_num_pairs (N) % birth_dim_1 % birth_value_1 % death_value_1 % ... % birth_dim_N % birth_value_N % death_value_N
                         (birth_values and death_values are of type double, the rest is int64_t) */
        template< typename Complex >
        void save_persistence_diagram( const std::string& filename,
                                       const inputs::abstract_weighted_cell_complex< Complex >& complex,
                                       const data_structures::distributed_vector< int64_t >& filtration_to_cell_map,
                                       const data_structures::write_once_column_array& reduced_columns,
                                       bool dualized = false,
                                       int64_t upper_dim = std::numeric_limits< int64_t >::max(),
                                       double upper_value = std::numeric_limits< double >::max(),
                                       bool without_top_dimension_essentials = true )
        {
            MPI_File file = mpi_utils::file_open_read_write( filename );

            const int64_t global_num_cols = complex.get_num_cells();
            const int64_t col_begin = element_distribution::get_local_begin( global_num_cols );
            const int64_t col_end = element_distribution::get_local_end( global_num_cols );

            // reserve space for preamble
            const int64_t preamble_length = 3 * sizeof( int64_t );
            MPI_Offset cur_file_size = preamble_length;
            MPI_File_set_size( file, cur_file_size );

            typedef std::pair< int64_t, std::pair< double, double > > diagram_entry_type;

            std::vector< diagram_entry_type > local_dims_and_pairs;
            std::vector< int64_t > birth_cell_queries;
            std::vector< int64_t > birth_cell_answers;
            std::vector< int64_t > death_cell_queries;
            std::vector< int64_t > death_cell_answers;
            std::vector< int64_t > birth_dim_queries;
            std::vector< int64_t > death_value_queries;
            std::vector< int64_t > birth_value_queries;
            std::vector< int64_t > birth_dim_answers;
            std::vector< double > birth_value_answers;
            std::vector< double > death_value_answers;
            std::vector< std::pair< int64_t, bool > > non_essential_cells;

            /// get maximum value of complex
            double max_value = std::numeric_limits< double >::max();
            if( upper_value == std::numeric_limits< double >::max() ) {
                max_value = complex.get_max_value();
            } else {
                double local_max_value_smaller_than = std::numeric_limits< double >::lowest( );
                for( int64_t idx = col_begin; idx < col_end; idx++ ) {
                    double value = complex.get_local_value( idx );
                    if( value <= upper_value )
                        local_max_value_smaller_than = value > local_max_value_smaller_than ? value : local_max_value_smaller_than;
                }
                std::vector< double > max_value_smaller_than_per_rank;
                mpi_utils::all_gather( local_max_value_smaller_than, max_value_smaller_than_per_rank );
                max_value = *std::max_element( max_value_smaller_than_per_rank.begin( ), max_value_smaller_than_per_rank.end( ) );
            }

            const int64_t max_dim = std::min( complex.get_max_dim( ), upper_dim );

            data_structures::distributed_vector< bool > is_cell_essential;
            is_cell_essential.init( global_num_cols );
            for( int64_t idx = col_begin; idx < col_end; idx++ )
                is_cell_essential.set_local_value( idx, true );

            const int64_t local_num_cells = col_end - col_begin;
            const int64_t num_chunks = ( local_num_cells + globals::DIPHA_BLOCK_SIZE ) / globals::DIPHA_BLOCK_SIZE;
            for( int64_t chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++ ) {
                const int64_t chunk_begin = col_begin + chunk_idx * globals::DIPHA_BLOCK_SIZE;
                const int64_t chunk_end = std::min( col_begin + ( chunk_idx + 1 ) * globals::DIPHA_BLOCK_SIZE, col_end );

                /// get birth cells indices
                birth_cell_queries.clear();
                for( int64_t idx = chunk_begin; idx < chunk_end; idx++ ) {
                    if( !reduced_columns.empty( idx ) ) {
                        if( dualized )
                            birth_cell_queries.push_back( reduced_columns.get_cell_index( idx ) );
                        else
                            birth_cell_queries.push_back( reduced_columns.front( idx ) );
                    }
                }
                birth_cell_answers.clear();
                filtration_to_cell_map.get_global_values( birth_cell_queries, birth_cell_answers );
                auto iterator_of_birth_cell_answers = birth_cell_answers.cbegin();

                /// get death cells indices
                death_cell_queries.clear();
                for( int64_t idx = chunk_begin; idx < chunk_end; idx++ ) {
                    if( !reduced_columns.empty( idx ) ) {
                        if( dualized )
                            death_cell_queries.push_back( reduced_columns.front( idx ) );
                        else
                            death_cell_queries.push_back( reduced_columns.get_cell_index( idx ) );
                    }
                }
                death_cell_answers.clear();
                filtration_to_cell_map.get_global_values( death_cell_queries, death_cell_answers );
                auto iterator_of_death_cell_answers = death_cell_answers.cbegin();

                ///  get birth dimension and birth / death values, and gather non_essential_cols
                birth_dim_queries.clear();
                death_value_queries.clear();
                birth_value_queries.clear();
                non_essential_cells.clear( );
                for( int64_t idx = chunk_begin; idx < chunk_end; idx++ ) {
                    if( !reduced_columns.empty( idx ) ) {
                        int64_t birth_cell = *iterator_of_birth_cell_answers++;
                        int64_t death_cell = *iterator_of_death_cell_answers++;
                        birth_dim_queries.push_back( birth_cell );
                        birth_value_queries.push_back( birth_cell );
                        death_value_queries.push_back( death_cell );
                        non_essential_cells.push_back( std::make_pair( birth_cell, false ) );
                        non_essential_cells.push_back( std::make_pair( death_cell, false ) );
                    }
                }
                complex.get_global_dims( birth_dim_queries, birth_dim_answers );
                auto iterator_of_birth_dim_answers = birth_dim_answers.cbegin();
                complex.get_global_values( birth_value_queries, birth_value_answers );
                auto iterator_of_birth_value_answers = birth_value_answers.cbegin();
                complex.get_global_values( death_value_queries, death_value_answers );
                auto iterator_of_death_value_answers = death_value_answers.cbegin();

                /// create local diagram entries
                for( int64_t idx = chunk_begin; idx < chunk_end; idx++ ) {
                    if( !reduced_columns.empty( idx ) ) {
                        int64_t birth_dim = *iterator_of_birth_dim_answers++;
                        double birth_value = *iterator_of_birth_value_answers++;
                        double death_value = *iterator_of_death_value_answers++;
                        if( birth_value != death_value && birth_value < max_value ) {
                            if( death_value <= max_value ) {
                                local_dims_and_pairs.push_back( std::make_pair( birth_dim, std::make_pair( birth_value, death_value ) ) );
                            } else {
                                int64_t shifted_dim = -birth_dim - 1;
                                local_dims_and_pairs.push_back( std::make_pair( shifted_dim, std::make_pair( birth_value, max_value ) ) );
                            }
                        }
                    }
                }

                /// update is_cell_essential
                is_cell_essential.set_global_values( non_essential_cells );
            }

            /// create essential triples
            for( int64_t idx = col_begin; idx < col_end; idx++ ) {
                if( is_cell_essential.get_local_value( idx ) == true ) {
                    int64_t dim = complex.get_local_dim( idx );
                    if( dim <= max_dim && ( !without_top_dimension_essentials || dim < max_dim ) ) {
                        int64_t shifted_dim = -dim - 1;
                        double value = complex.get_local_value( idx );
                        if( value <= max_value )
                            local_dims_and_pairs.push_back( std::make_pair( shifted_dim, std::make_pair( value, max_value ) ) );
                    }
                }
            }

            int64_t local_num_pairs = local_dims_and_pairs.size( );
            std::vector< int64_t > local_num_pairs_per_rank;
            mpi_utils::all_gather( local_num_pairs, local_num_pairs_per_rank );
            std::vector< int64_t > cum_sum_local_num_pairs( mpi_utils::get_num_processes( ) + 1 );
            cum_sum_local_num_pairs.front( ) = 0;
            std::partial_sum( local_num_pairs_per_rank.begin( ), local_num_pairs_per_rank.end( ), cum_sum_local_num_pairs.begin( ) + 1 );
            int64_t global_num_pairs = cum_sum_local_num_pairs.back( );

            const int64_t entry_size = sizeof( diagram_entry_type );
            MPI_Offset file_offset = cur_file_size + cum_sum_local_num_pairs[ mpi_utils::get_rank( ) ] * entry_size;

            /// reserve space for pairs
            cur_file_size = preamble_length + global_num_pairs * entry_size;
            MPI_File_set_size( file, cur_file_size );

            /// now write pairs to file
            mpi_utils::file_write_at_vector( file, file_offset, local_dims_and_pairs );

            /// write preamble to file
            if( mpi_utils::is_root() ) {
                mpi_utils::file_write_at( file, 0, file_types::DIPHA );
                mpi_utils::file_write_at( file, sizeof( int64_t ), file_types::PERSISTENCE_DIAGRAM );
                mpi_utils::file_write_at( file, 2 * sizeof( int64_t ), global_num_pairs );
            }

            /// now sort the pairs to simplify regression tests etc.
            MPI_Barrier( MPI_COMM_WORLD );

            /// be careful with tiny files: psort needs all processes to have work to do ...
            if( global_num_pairs < mpi_utils::get_num_processes( ) * mpi_utils::get_num_processes( ) ) {
                if( mpi_utils::is_root( ) ) {
                    std::vector< diagram_entry_type > global_dims_and_pairs;
                    MPI_Offset offset = 3 * sizeof( int64_t );
                    mpi_utils::file_read_at_vector( file, offset, global_num_pairs, global_dims_and_pairs );
                    std::sort( global_dims_and_pairs.begin( ), global_dims_and_pairs.end( ) );
                    mpi_utils::file_write_at_vector( file, offset, global_dims_and_pairs );
                }
            } else {
                const int64_t local_begin = element_distribution::get_local_begin( global_num_pairs );
                const int64_t local_num_pairs = element_distribution::get_local_end( global_num_pairs ) - local_begin;

                std::vector< diagram_entry_type > local_dims_and_pairs;

                MPI_Offset offset = sizeof(int64_t)* 3 + local_begin * sizeof( diagram_entry_type );
                mpi_utils::file_read_at_vector( file, offset, local_num_pairs, local_dims_and_pairs );

                // psort unfortunately uses long for the size. This will cause problems on Win64 for large data
                std::vector< long > cell_distribution;
                int num_processes = mpi_utils::get_num_processes( );
                for( int cur_rank = 0; cur_rank < num_processes; cur_rank++ )
                    cell_distribution.push_back( (long)( element_distribution::get_local_end( global_num_pairs, cur_rank ) - element_distribution::get_local_begin( global_num_pairs, cur_rank ) ) );
                p_sort::parallel_sort( local_dims_and_pairs.begin( ), local_dims_and_pairs.end( ), cell_distribution.data( ), MPI_COMM_WORLD );

                // need to make sure that above read has completeted before we overwrite
                MPI_Barrier( MPI_COMM_WORLD );
                mpi_utils::file_write_at_vector( file, offset, local_dims_and_pairs );
            }

            MPI_File_close( &file );
        }
    }
}