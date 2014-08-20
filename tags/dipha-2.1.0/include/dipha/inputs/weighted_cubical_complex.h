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
    namespace inputs {
        class weighted_cubical_complex : public abstract_weighted_cell_complex< weighted_cubical_complex > {

            friend abstract_weighted_cell_complex< weighted_cubical_complex >;

            // defining state of the object
        protected:
            std::vector< int64_t > lattice_resolution;
            std::vector< double > vertex_values;

            // derived quantities
        protected:
            std::vector< int64_t > cum_resolution_product;
            std::vector< int64_t > cubical_complex_resolution;
            int64_t vertex_values_begin;
            int64_t vertex_values_end;

            // implementation of abstract_weighted_cell_complex interface
        protected:
            int64_t _get_num_cells() const
            {
                int64_t num_cells = 1;
                for( int64_t cur_dim = 0; cur_dim < get_max_dim(); cur_dim++ )
                    num_cells *= cubical_complex_resolution[ cur_dim ];
                return num_cells;
            }

            int64_t _get_max_dim() const { return lattice_resolution.size(); }

            int64_t _get_local_dim( int64_t idx ) const
            {
                static std::vector< int64_t > temp_tuple;
                index_to_tuple( idx, temp_tuple );
                int64_t cell_dim = 0;
                for( int64_t cur_dim = 0; cur_dim < get_max_dim(); cur_dim++ )
                    cell_dim += temp_tuple[ cur_dim ] % 2;
                return cell_dim;
            }

            double _get_local_value( int64_t idx ) const
            {
                static std::vector< int64_t > temp_tuple;
                index_to_tuple( idx, temp_tuple );

                static std::vector< int64_t > target_dimensions;
                target_dimensions.clear();
                for( int64_t cur_dim = 0; cur_dim < get_max_dim(); cur_dim++ )
                if( temp_tuple[ cur_dim ] % 2 )
                    target_dimensions.push_back( cur_dim );

                static std::vector< int64_t > vertex_tuple;
                vertex_tuple.clear();
                double max_vertex_value = std::numeric_limits< double >::lowest();
                int64_t num_vertices = 2 << target_dimensions.size();
                for( int64_t k = 0; k < num_vertices; k++ ) {
                    vertex_tuple = temp_tuple;
                    for( int64_t j = 0; j < (int64_t)target_dimensions.size(); j++ )
                        vertex_tuple[ target_dimensions[ j ] ] += -1 + 2 * ( k / ( 2 << j ) % 2 == 0 );
                    double vertex_value = vertex_values[ vertex_tuple_to_lattice_index( vertex_tuple ) - vertex_values_begin ];
                    max_vertex_value = vertex_value > max_vertex_value ? vertex_value : max_vertex_value;
                }

                return max_vertex_value;
            }

            void _get_local_boundary( int64_t idx, std::vector< int64_t >& boundary ) const
            {
                boundary.clear();
                static std::vector< int64_t > temp_tuple;
                index_to_tuple( idx, temp_tuple );
                for( int64_t cur_dim = 0; cur_dim < get_max_dim(); cur_dim++ ) {
                    if( temp_tuple[ cur_dim ] % 2 ) {
                        temp_tuple[ cur_dim ]--;
                        boundary.push_back( tuple_to_index( temp_tuple ) );
                        temp_tuple[ cur_dim ] += 2;
                        boundary.push_back( tuple_to_index( temp_tuple ) );
                        temp_tuple[ cur_dim ]--;
                    }
                }
            }

            void _get_local_coboundary( int64_t idx, std::vector< int64_t >& coboundary ) const
            {
                coboundary.clear();
                static std::vector< int64_t > temp_tuple;
                index_to_tuple( idx, temp_tuple );
                for( int64_t cur_dim = 0; cur_dim < get_max_dim(); cur_dim++ ) {
                    if( temp_tuple[ cur_dim ] % 2 == 0 ) {
                        temp_tuple[ cur_dim ]--;
                        if( temp_tuple[ cur_dim ] >= 0 ) {
                            coboundary.push_back( tuple_to_index( temp_tuple ) );
                        }
                        temp_tuple[ cur_dim ] += 2;
                        if( temp_tuple[ cur_dim ] < cubical_complex_resolution[ cur_dim ] ) {
                            coboundary.push_back( tuple_to_index( temp_tuple ) );
                        }

                        temp_tuple[ cur_dim ]--;
                    }
                }
            }

            // Loads the weighted_cubical_complex from given file in binary format -- all symbols are 64 bit wide
            // Format: file_types::DIPHA % file_types::IMAGE_DATA % num_vertices % max_dim % lattice_resolution_1 % ... % lattice_resolution_dim % value_1 % ... % value_num_vertices
            void _load_binary( MPI_File file,
                               int64_t upper_dim = std::numeric_limits< int64_t >::max( ),
                               double upper_value = std::numeric_limits< double >::max( ) )
            {
                // read preamble
                std::vector< int64_t > preamble;
                mpi_utils::file_read_at_vector( file, 0, 4, preamble );

                int64_t global_num_vertices = preamble[ 2 ];
                int64_t max_dim = preamble[ 3 ];

                MPI_Offset lattice_resolution_offset = preamble.size() * sizeof( int64_t );
                mpi_utils::file_read_at_vector( file, lattice_resolution_offset, max_dim, lattice_resolution );

                cubical_complex_resolution.resize( lattice_resolution.size( ) );
                for( int64_t cur_dim = 0; cur_dim < get_max_dim( ); cur_dim++ )
                    cubical_complex_resolution[ cur_dim ] = 2 * lattice_resolution[ cur_dim ] - 1;

                cum_resolution_product.resize( get_max_dim( ) - 1 );
                cum_resolution_product[ 0 ] = cubical_complex_resolution[ 0 ];
                for( int64_t cur_dim = 1; cur_dim < get_max_dim( ) - 1; cur_dim++ )
                    cum_resolution_product[ cur_dim ] = cum_resolution_product[ cur_dim - 1 ] * cubical_complex_resolution[ cur_dim ];

                vertex_values_begin = std::numeric_limits< int64_t >::max();
                vertex_values_end = -1;
                std::vector< int64_t > target_dimensions;
                std::vector< int64_t > temp_tuple;
                std::vector< int64_t > vertex_tuple;
                const int64_t local_begin = element_distribution::get_local_begin( get_num_cells() );
                const int64_t local_end = element_distribution::get_local_end( get_num_cells() );
                for( int64_t local_idx = local_begin; local_idx < local_end; local_idx++ ) {
                    index_to_tuple( local_idx, temp_tuple );
                    target_dimensions.clear();
                    for( int64_t cur_dim = 0; cur_dim < get_max_dim(); cur_dim++ ) {
                        if( temp_tuple[ cur_dim ] % 2 )
                            target_dimensions.push_back( cur_dim );
                    }

                    vertex_tuple.clear();
                    int64_t num_vertices = 2 << target_dimensions.size();
                    for( int64_t k = 0; k < num_vertices; k++ ) {
                        vertex_tuple = temp_tuple;
                        for( int64_t j = 0; j < (int64_t)target_dimensions.size(); j++ )
                            vertex_tuple[ target_dimensions[ j ] ] += -1 + 2 * ( k / ( 2 << j ) % 2 == 0 );
                        int64_t lattice_index = vertex_tuple_to_lattice_index( vertex_tuple );
                        vertex_values_begin = vertex_values_begin < lattice_index ? vertex_values_begin : lattice_index;
                        vertex_values_end = vertex_values_end > lattice_index ? vertex_values_end : lattice_index;
                    }
                }
                vertex_values_end += 1;

                int64_t local_num_vertices = vertex_values_end - vertex_values_begin;
                MPI_Offset vertex_values_offset = ( preamble.size() + max_dim + vertex_values_begin ) * sizeof( int64_t );
                mpi_utils::file_read_at_vector( file, vertex_values_offset, local_num_vertices, vertex_values );
            }

            // overide default implemantations in abstract_weighted_cell_complex to improve performance
        protected:
            double _get_max_value( ) const
            {
                const int64_t local_begin = element_distribution::get_local_begin( get_num_cells( ) );
                const int64_t local_end = element_distribution::get_local_end( get_num_cells( ) );
                double local_max_value = *std::max_element( vertex_values.begin( ), vertex_values.end( ) );
                std::vector< double > max_value_per_rank;
                mpi_utils::all_gather( local_max_value, max_value_per_rank );
                return *std::max_element( max_value_per_rank.begin( ), max_value_per_rank.end( ) );
            }

            void _get_global_dims( const std::vector< int64_t >& queries,
                                   std::vector< int64_t >& answers ) const
            {
                answers.resize( queries.size() );
                for( int64_t idx = 0; idx < (int64_t)queries.size(); idx++ )
                    answers[ idx ] = _get_local_dim( queries[ idx ] );
            }

            void _get_global_boundaries( const std::vector< int64_t >& queries,
                                         data_structures::write_once_array_of_arrays< int64_t >& answers ) const
            {
                std::vector< int64_t > boundary;
                for( int64_t idx = 0; idx < (int64_t)queries.size(); idx++ ) {
                    _get_local_boundary( queries[ idx ], boundary );
                    answers.set( idx, boundary.begin(), boundary.end() );
                }
            }

            void _get_global_coboundaries( const std::vector< int64_t >& queries,
                                           data_structures::write_once_array_of_arrays< int64_t >& answers ) const
            {
                std::vector< int64_t > coboundary;
                for( int64_t idx = 0; idx < (int64_t)queries.size(); idx++ ) {
                    _get_local_coboundary( queries[ idx ], coboundary );
                    answers.set( idx, coboundary.begin(), coboundary.end() );
                }
            }

            // internal helper methods
        protected:
            int64_t tuple_to_index( const std::vector< int64_t >& tuple ) const
            {
                int64_t index = tuple[ get_max_dim() - 1 ];
                for( int64_t cur_dim = get_max_dim() - 2; cur_dim >= 0; cur_dim-- ) {
                    index *= cubical_complex_resolution[ cur_dim ];
                    index += tuple[ cur_dim ];
                }
                return index;
            }

            int64_t vertex_tuple_to_lattice_index( const std::vector< int64_t >& tuple ) const
            {
                int64_t index = tuple[ get_max_dim() - 1 ] / 2;
                for( int64_t cur_dim = get_max_dim() - 2; cur_dim >= 0; cur_dim-- ) {
                    index *= lattice_resolution[ cur_dim ];
                    index += tuple[ cur_dim ] / 2;
                }
                return index;
            }

            void index_to_tuple( int64_t idx, std::vector< int64_t >& tuple ) const
            {
                tuple.resize( get_max_dim() );
                for( int64_t cur_dim = get_max_dim() - 1; cur_dim > 0; cur_dim-- ) {
                    tuple[ cur_dim ] = idx / cum_resolution_product[ cur_dim - 1 ];
                    idx -= tuple[ cur_dim ] * cum_resolution_product[ cur_dim - 1 ];
                }
                tuple[ 0 ] = idx;
            }
        };
    }
}
