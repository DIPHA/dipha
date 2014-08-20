/*  Copyright 2014 IST Austria

Contributed by: Michael Kerber

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
        class sparse_rips_complex : public abstract_weighted_cell_complex< sparse_rips_complex > {

            friend abstract_weighted_cell_complex< sparse_rips_complex >;

            // defining state of the object
        protected:
            int64_t _m_upper_dim;
            int64_t _m_no_points;

            std::vector<std::vector<double> > _m_distance_matrix;

            double _m_threshold;

            // derived quantities
        protected:

            std::vector< int64_t > _m_full_indices_in_range;

            std::vector<std::vector<int64_t> > _m_binomials;

            std::vector< int64_t > _m_breakpoints_local_indices;

            int64_t _m_num_elements;

            // contains the indices where a new dimension starts
            std::vector<int64_t> _m_breakpoints;

        protected:
            // Loads the complete_rips_complex from given file in binary format -- all symbols are 64 bit wide
            // Format: magic_number % file_type % num_points % max_dim of complex % values of distance matrix
            // TODO: Complete matrix, or just upper part? For now, I did the whole matrix
            void _load_binary( MPI_File file, 
                               int64_t upper_dim = std::numeric_limits< int64_t >::max(),
                               double upper_value = std::numeric_limits< double >::max( ) )
            {
                // read preamble
                std::vector< int64_t > preamble;
                mpi_utils::file_read_at_vector( file, 0, 3, preamble );

                _m_no_points = preamble[ 2 ];
                
                if( upper_dim == std::numeric_limits< int64_t >::max( ) ) {
                    mpi_utils::error_printer_if_root( ) << "No upper_dim specified for distance_matrix data!";
                    MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
                } else {
                    _m_upper_dim = upper_dim;
                }

                _m_threshold = 2 * upper_value;

                int64_t matrix_size = _m_no_points * _m_no_points;

                std::vector< double > distances;
                MPI_Offset offset = preamble.size( ) * sizeof(int64_t);
                mpi_utils::file_read_at_vector( file, offset, matrix_size, distances );


                _m_distance_matrix.resize( _m_no_points );
                for( int64_t i = 0; i < _m_no_points; i++ ) {
                    _m_distance_matrix[ i ].resize( _m_no_points );
                }

                for( int64_t i = 0; i < matrix_size; i++ ) {
                    _m_distance_matrix[ i / _m_no_points ][ i % _m_no_points ] = distances[ i ];
                }

                _precompute_binomials();

                _precompute_breakpoints();

                _compute_sparse_indices();

            }

        protected:

            // Local version, in order to render hashing possible, possibly
            int get_process_id_for_full_index( int64_t full_idx )
            {
                return dipha::element_distribution::get_rank( num_full_indices(), full_idx );
            }


            int64_t num_full_indices() const
            {
                return _m_breakpoints[ _m_upper_dim + 1 ];
            }

            // Note: This currently scales very badly with the number of points
            // (TODO: Sparse list representation?)
            // scan_full_column=false means: only look for the coboundaries such that
            // the new vertex has the largest index among all indices
            void _get_local_coboundary_full_index( int64_t full_idx,
                                                   std::vector<int64_t>& coboundary,
                                                   bool scan_full_column = false ) const
            {

	      //std::cout << "compute_coboundary for " << full_idx << " " << std::endl;

                coboundary.clear();

		int64_t k = _get_local_dim_full_index(full_idx);

		if(k==_m_upper_dim) {
		  return;
		}

                std::vector<int64_t> vertex_indices;
                conversion( full_idx, std::back_inserter( vertex_indices ) );

                int64_t start_idx_for_scan = scan_full_column ? 0 : vertex_indices.front();

                int64_t rows_to_check = vertex_indices.size();

                for( int64_t coboundary_candidate = start_idx_for_scan; coboundary_candidate < _m_no_points; coboundary_candidate++ ) {

                    bool makes_a_coboundary = true;

                    for( int row_idx = 0; row_idx < rows_to_check; row_idx++ ) {
                        int64_t row = vertex_indices[ row_idx ];
                        if( row == coboundary_candidate || _m_distance_matrix[ row ][ coboundary_candidate ] > _m_threshold ) {
                            makes_a_coboundary = false;
                            //std::cout << "Column " << coboundary_candidate << "is not a coboundary for " << full_idx << std::endl;
                            break;
                        }
                    }

                    if( makes_a_coboundary ) {

                        int64_t full_idx_of_new_coboundary;
                        conversion_with_extra_index( vertex_indices.begin(), vertex_indices.end(), coboundary_candidate, full_idx_of_new_coboundary );
                        coboundary.push_back( full_idx_of_new_coboundary );
                    }

                }
            }



            void _compute_sparse_indices()
            {

	      //std::cout << "Compute sparse_indices.." << std::flush;

                std::vector<std::vector<int64_t> > full_indices_to_handle_per_dimension;
                full_indices_to_handle_per_dimension.resize( _m_upper_dim + 1 );


                int num_processes = dipha::mpi_utils::get_num_processes();

                // Assign the vertices to the right processes (normally, the 0th process should get all of them without hashing)
                std::vector<int64_t> simplex_of_dim_0;
                simplex_of_dim_0.push_back( -1 );

                int64_t full_index_of_vertex;

                for( int point_idx = 0; point_idx < _m_no_points; point_idx++ ) {
                    simplex_of_dim_0[ 0 ] = point_idx;
                    conversion_with_skip( simplex_of_dim_0.begin(), simplex_of_dim_0.end(), simplex_of_dim_0.end(), full_index_of_vertex );

                    int64_t process_id_for_vertex = get_process_id_for_full_index( full_index_of_vertex );

                    if( process_id_for_vertex == dipha::mpi_utils::get_rank() ) {
                        full_indices_to_handle_per_dimension[ 0 ].push_back( full_index_of_vertex );
                    }
                }

                std::copy( full_indices_to_handle_per_dimension[ 0 ].begin(), full_indices_to_handle_per_dimension[ 0 ].end(),
                           std::back_inserter( _m_full_indices_in_range ) );

                // Go in increasing dimension...
                for( int dim = 0; dim < _m_upper_dim; dim++ ) {

                    // find out the simplices in codimension one
                    std::vector<int64_t> full_index_of_codimension_one_buffer;

                    //std::cout << "NOW " << local_process_idx << " " << _m_coboundaries.size() << std::endl;

                    for( auto full_index_it = full_indices_to_handle_per_dimension[ dim ].begin();
                         full_index_it != full_indices_to_handle_per_dimension[ dim ].end();
                         full_index_it++ ) {
                        std::vector<int64_t> coboundary;
                        _get_local_coboundary_full_index( *full_index_it, coboundary, false );
                        std::copy( coboundary.begin(), coboundary.end(), std::back_inserter( full_index_of_codimension_one_buffer ) );
                    }

                    std::vector<std::vector<int64_t> > codimension_one_full_indices_send_buffer;
                    codimension_one_full_indices_send_buffer.resize( num_processes );

                    // Send the newly found indices to the appropriate indices

                    MPI_Barrier( MPI_COMM_WORLD );

                    for( auto new_idx_it = full_index_of_codimension_one_buffer.begin();
                         new_idx_it != full_index_of_codimension_one_buffer.end();
                         new_idx_it++ ) {

                        int responsible_process = get_process_id_for_full_index( *new_idx_it );
                        codimension_one_full_indices_send_buffer[ responsible_process ].push_back( *new_idx_it );
                    }
                    std::vector< MPI_Request > queries_requests;
                    for( int process_id = 0; process_id < num_processes; process_id++ ) {
                        mpi_utils::non_blocking_send_vector( codimension_one_full_indices_send_buffer[ process_id ], process_id, mpi_utils::MSG_REPORT_INDICES_OF_SPARSE_RIPS, queries_requests );
                    }


                    // ...and receive the results from other processes to handle the next dimension

                    std::vector<std::vector<int64_t> > codimension_one_full_indices_receive_buffer;
                    codimension_one_full_indices_receive_buffer.resize( num_processes );

                    for( int process_id = 0; process_id < num_processes; process_id++ ) {
                        int source = ( dipha::mpi_utils::get_rank() + process_id ) % num_processes;
                        mpi_utils::receive_vector( codimension_one_full_indices_receive_buffer[ source ], source, mpi_utils::MSG_REPORT_INDICES_OF_SPARSE_RIPS );
                        std::copy( codimension_one_full_indices_receive_buffer[ source ].begin(),
                                   codimension_one_full_indices_receive_buffer[ source ].end(),
                                   std::back_inserter( full_indices_to_handle_per_dimension[ dim + 1 ] ) );
                    }

                    std::sort( full_indices_to_handle_per_dimension[ dim + 1 ].begin(), full_indices_to_handle_per_dimension[ dim + 1 ].end() );

                    std::copy( full_indices_to_handle_per_dimension[ dim + 1 ].begin(), full_indices_to_handle_per_dimension[ dim + 1 ].end(),
                               std::back_inserter( _m_full_indices_in_range ) );

                    MPI_Waitall( (int)queries_requests.size(), queries_requests.data(), MPI_STATUSES_IGNORE );

                }

                // Now, tell everyone how many simplices are you responsible for:
                int64_t no_of_local_simplices = _m_full_indices_in_range.size();

                // Silly: Create vector of size one to use send_vector from mpi_utils
                std::vector<int64_t> wrap_no_of_local_simplices;

                wrap_no_of_local_simplices.push_back( no_of_local_simplices );

                std::vector< MPI_Request > queries_requests;

                for( int process_id = 0; process_id < num_processes; process_id++ ) {
                    mpi_utils::non_blocking_send_vector( wrap_no_of_local_simplices, process_id, mpi_utils::MSG_REPORT_INDICES_OF_SPARSE_RIPS, queries_requests );
                }

                // And gather the information of each process
                _m_breakpoints_local_indices.push_back( 0 );
                std::vector<int64_t> no_of_local_simplices_from_process;
                for( int process_id = 0; process_id < num_processes; process_id++ ) {
                    // Silly again: Receive using a vector of size one
                    mpi_utils::receive_vector( no_of_local_simplices_from_process, process_id, mpi_utils::MSG_REPORT_INDICES_OF_SPARSE_RIPS );
                    assert( no_of_local_simplices_from_process.size() == 1 );
                    _m_breakpoints_local_indices.push_back( no_of_local_simplices_from_process[ 0 ] );
                }
                // Finally, sum up the simplices
                for( int64_t process_id = 1; process_id <= num_processes; process_id++ ) {
                    _m_breakpoints_local_indices[ process_id ] += _m_breakpoints_local_indices[ process_id - 1 ];
                }
                _m_num_elements = _m_breakpoints_local_indices.back();

                //std::cout << "Found " << _m_num_elements << " simplices" << std::endl;

                MPI_Waitall( (int)queries_requests.size(), queries_requests.data(), MPI_STATUSES_IGNORE );

                // Now, each process knows that total number of simplices, and therefore, the range of repsonsibility for each process

                // Send the full indices to the responsible process
                std::vector< std::vector< int64_t > > send_buffers;
                send_buffers.resize( num_processes );

                for( int64_t local_sparse_idx = 0; local_sparse_idx < ( int64_t)_m_full_indices_in_range.size( ); local_sparse_idx++ )
                {
                    int64_t sparse_idx = _m_breakpoints_local_indices[ dipha::mpi_utils::get_rank() ] + local_sparse_idx;

                    int process_responsible_for_current = dipha::element_distribution::get_rank( _m_num_elements, sparse_idx );
                    //std::cout << dipha::mpi_utils::get_rank() << "Responsible is " << _m_num_elements << " " << local_sparse_idx << " " << sparse_idx << " " <<  process_responsible_for_current << std::endl;
                    send_buffers[ process_responsible_for_current ].push_back( _m_full_indices_in_range[ local_sparse_idx ] );
                }

                queries_requests.clear();
                for( int target = 0; target < num_processes; target++ )
                    dipha::mpi_utils::non_blocking_send_vector( send_buffers[ target ], target, dipha::mpi_utils::MSG_SPARSE_RIPS_DISTRIBUTE_SIMPLICES, queries_requests );

                std::vector< std::vector< int64_t> > queries_buffer;
                queries_buffer.resize( num_processes );
                for( int idx = 0; idx < num_processes; idx++ ) {
                    int source = ( dipha::mpi_utils::get_rank() + idx ) % num_processes;
                    dipha::mpi_utils::receive_vector( queries_buffer[ source ], source, dipha::mpi_utils::MSG_SPARSE_RIPS_DISTRIBUTE_SIMPLICES );
                }

                MPI_Waitall( (int)queries_requests.size(), queries_requests.data(), MPI_STATUSES_IGNORE );

                // queries buffer contains the indices of the process, store them (overwrite the _m_full_indices_in_range)
                _m_full_indices_in_range.clear();
                for( int64_t proc_idx = 0; proc_idx < ( int64_t)queries_buffer.size( ); proc_idx++ )
                {
                    std::copy( queries_buffer[ proc_idx ].begin(), queries_buffer[ proc_idx ].end(), std::back_inserter( _m_full_indices_in_range ) );
                }
                // sort the indices
                std::sort( _m_full_indices_in_range.begin(), _m_full_indices_in_range.end() );

                /*
                        for( int64_t i = 0; i < _m_full_indices_in_range.size(); i++ ) {
                        std::cout << dipha::mpi_utils::get_rank() << " " << _m_full_indices_in_range[ i ] << std::endl;
                        }
                        */

                if( _m_full_indices_in_range.empty() ) {
                    std::cout << std::endl << "Error: Complex too small for number of processes used!" << std::endl;
                    MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
                }

                // Finally, tell everyone the maximal global index in the range (so that everyone can ask at the right place)
                int64_t maximal_global_idx = _m_full_indices_in_range.back();
                // Silly: Create vector of size one to use send_vector from mpi_utils
                std::vector<int64_t> wrap_maximal_global_idx;

                wrap_maximal_global_idx.push_back( maximal_global_idx );

                queries_requests.clear();

                for( int process_id = 0; process_id < num_processes; process_id++ ) {
                    mpi_utils::non_blocking_send_vector( wrap_maximal_global_idx, process_id, dipha::mpi_utils::MSG_REPORT_MAXIMAL_GLOBAL_INDEX_IN_SPARSE_RIPS, queries_requests );
                }

                // And gather the information of each process
                // Overwrite _m_breakpoints_local_indices
                _m_breakpoints_local_indices.clear();
                _m_breakpoints_local_indices.push_back( 0 );
                std::vector<int64_t> breakpoint_from_process;
                for( int process_id = 0; process_id < num_processes; process_id++ ) {
                    // Silly again: Receive using a vector of size one
                    dipha::mpi_utils::receive_vector( breakpoint_from_process, process_id, dipha::mpi_utils::MSG_REPORT_MAXIMAL_GLOBAL_INDEX_IN_SPARSE_RIPS );
                    assert( breakpoint_from_process.size() == 1 );
                    _m_breakpoints_local_indices.push_back( 1 + breakpoint_from_process[ 0 ] );
                }
                /*
                        std::cout << "Process " << dipha::mpi_utils::get_rank() << " has " << _m_full_indices_in_range.size() << " simplices" << std::endl;
                        for( int64_t i = 0; i <= num_processes; i++ ) {
                        std::cout << dipha::mpi_utils::get_rank() << " breakpoint " << i << " " << _m_breakpoints_local_indices[ i ] << std::endl;
                        }
                        */

		//std::cout << "done" << std::endl;

            }

            int64_t _get_max_dim() const
            {
                return _m_upper_dim;
            }


            int64_t _get_num_cells() const
            {
                return _m_num_elements;
            }


            int64_t get_locally_sparse_from_full_index( int64_t full_idx ) const
            {
                auto it = std::lower_bound( _m_full_indices_in_range.begin(), _m_full_indices_in_range.end(), full_idx );
                assert( it != _m_full_indices_in_range.end() );
                assert( *it == full_idx );
                int64_t result = dipha::element_distribution::get_local_begin( _m_num_elements, dipha::mpi_utils::get_rank() ) + ( it - _m_full_indices_in_range.begin() );
                //std::cout << "full to sparse: " << dipha::mpi_utils::get_rank() << " " << full_idx << " " << result << std::endl;
                return result;
            }

            int64_t get_locally_full_from_sparse_index( int64_t sparse_idx ) const
            {

                int process_id = dipha::mpi_utils::get_rank();

                int64_t local_begin = dipha::element_distribution::get_local_begin( _m_num_elements, process_id );
                /*
                        if( sparse_idx >= dipha::element_distribution::get_local_end( _m_num_elements, process_id ) ) {
                        std::cout << mpi_utils::get_rank() << " " << "Got sparse_idx " << sparse_idx << ", but my range is " << local_begin << ", " << dipha::element_distribution::get_local_end( _m_num_elements, process_id ) << std::endl;
                        }
                        */
                assert( sparse_idx >= local_begin );
                assert( sparse_idx < dipha::element_distribution::get_local_end( _m_num_elements, process_id ) );

                sparse_idx -= local_begin;

                return _m_full_indices_in_range[ sparse_idx ];
            }

            int64_t _get_local_dim( int64_t idx ) const
            {
                int64_t full_idx = get_locally_full_from_sparse_index( idx );
                return _get_local_dim_full_index( full_idx );
            }

            int64_t _get_local_dim_full_index( int64_t idx ) const
            {
                int64_t k = 0;
                while( k < _m_upper_dim + 1 && _m_breakpoints[ k + 1 ] <= idx ) {
                    k++;
                }
                return k;
            }

            double _get_local_value( int64_t idx ) const
            {
                int64_t full_idx = get_locally_full_from_sparse_index( idx );
                return _get_local_value_full_index( full_idx );
            }


            double _get_local_value_full_index( int64_t idx ) const
            {
                static std::vector<int64_t> points;
                points.clear();
                conversion( idx, std::back_inserter( points ) );
                return diameter( points.begin(), points.end() );
            }

            /* Not possible to determine locally (at least in this way)
                void _get_local_boundary( int64_t idx, std::vector< int64_t >& boundary ) const
                {
                int64_t full_idx = get_locally_full_from_sparse_index(idx);
                return _get_local_boundary_full_index( full_idx, boundary );
                }
                */

            void _get_global_boundaries( const std::vector< int64_t >& queries,
                                         data_structures::write_once_array_of_arrays< int64_t >& answers ) const
            {
                return _get_global_co_boundaries( queries, answers, false );
            }

            void _get_global_coboundaries( const std::vector< int64_t >& queries,
                                           data_structures::write_once_array_of_arrays< int64_t >& answers ) const
            {
                return _get_global_co_boundaries( queries, answers, true );
            }

            int process_for_full_index( int64_t full_idx ) const
            {
                auto it = std::upper_bound( _m_breakpoints_local_indices.begin(), _m_breakpoints_local_indices.end(), full_idx );
                //std::cout << "pffi returns " << full_idx << " -> " << ( it - _m_breakpoints_local_indices.begin() - 1 ) << std::endl;
                return (int)(it - _m_breakpoints_local_indices.begin() - 1);
            }

            void _get_global_co_boundaries( const std::vector< int64_t >& queries,
                                            data_structures::write_once_array_of_arrays< int64_t >& answers,
                                            bool dual ) const
            {
                int num_processes = dipha::mpi_utils::get_num_processes();

                std::vector< std::vector< int64_t > > queries_buffer;
                element_distribution::scatter_queries( queries, _m_num_elements, queries_buffer );

                // process queries (This is a two step approach)

                // first, we compute which full indices are needed at all
#if 1
                std::unordered_map< int64_t, int64_t > full_to_sparse_indices_in_co_boundary;
#else
		std::map< int64_t, int64_t > full_to_sparse_indices_in_co_boundary;
#endif
                std::vector< int64_t > boundary;
                for( int source = 0; source < num_processes; source++ ) {
                    for( int64_t idx = 0; idx < (int64_t)queries_buffer[ source ].size(); idx++ ) {
                        //std::cout << "New el_idx (sparse): " << dipha::mpi_utils::get_rank() << " " << queries_buffer[ source ][ idx ] << std::endl;
                        int64_t full_idx = get_locally_full_from_sparse_index( queries_buffer[ source ][ idx ] );
                        if( dual ) {

                            //std::cout << "YYY: Full Index: " << full_idx << std::endl;
                            _get_local_coboundary_full_index( full_idx, boundary, true );
                            //std::cout << "YYY: Result: ";
                            /*
                            for(int i=0;i<boundary.size();i++) {
                            std::cout << boundary[i] << " ";
                            }
                            std::cout << std::endl;
                            */
                        } else {
                            _get_local_boundary_full_index( full_idx, boundary );
                        }
                        for( auto it = boundary.begin(); it != boundary.end(); it++ ) {
                            full_to_sparse_indices_in_co_boundary[ *it ] = -1; // dummy entry for now
                        }
                    }
                }

                // Next, we need to ask for the sparse indices of all collected full indices 
                std::vector< std::vector< int64_t > > inner_queries_send_buffer;
                std::vector< std::vector< int64_t > > inner_queries_recv_buffer;
                // std::vector< std::vector< int64_t > > inner_answers_send_buffer;
                std::vector< std::vector< int64_t > > inner_answers_recv_buffer;

                inner_queries_send_buffer.resize( num_processes );

                //std::cout << "KEYS " << std::flush;

                for( auto map_it = full_to_sparse_indices_in_co_boundary.begin();
                     map_it != full_to_sparse_indices_in_co_boundary.end();
                     map_it++ ) {
                    int64_t key = map_it->first;
                    //std::cout << key << " " << std::flush;
                    assert( map_it->second == -1 );
                    int process_for_current = process_for_full_index( key );
                    //std::cout << "To process " << process_for_current << "! " << std::flush;
                    inner_queries_send_buffer[ process_for_current ].push_back( key );
                }
                //std::cout << std::endl;

                std::vector< MPI_Request > queries_requests;

                // Send the queries now to ask for the sparse indices
                for( int process_id = 0; process_id < num_processes; process_id++ ) {
                    mpi_utils::non_blocking_send_vector( inner_queries_send_buffer[ process_id ], process_id, mpi_utils::MSG_QUERY_SPARSE_INDICES, queries_requests );
                }

                // Receive the query results 

                inner_queries_recv_buffer.resize( num_processes );

                for( int process_id = 0; process_id < num_processes; process_id++ ) {
                    int source = ( dipha::mpi_utils::get_rank( ) + process_id ) % num_processes;
                    mpi_utils::receive_vector( inner_queries_recv_buffer[ source ], source, mpi_utils::MSG_QUERY_SPARSE_INDICES );
                }

                MPI_Waitall( (int)queries_requests.size(), queries_requests.data(), MPI_STATUSES_IGNORE );


                // Sent out the answers
                for( int process_id = 0; process_id < num_processes; process_id++ ) {
                    int source = ( dipha::mpi_utils::get_rank() + process_id ) % num_processes;
                    for( int64_t query_buf_idx = 0; query_buf_idx < (int64_t)inner_queries_recv_buffer[ source ].size( ); query_buf_idx++ )
                    {
                        inner_queries_recv_buffer[ source ][ query_buf_idx ] = get_locally_sparse_from_full_index( inner_queries_recv_buffer[ source ][ query_buf_idx ] );
                    }
                    mpi_utils::non_blocking_send_vector( inner_queries_recv_buffer[ source ], source, mpi_utils::MSG_QUERY_SPARSE_INDICES, queries_requests );
                }

                inner_answers_recv_buffer.resize( num_processes );
                // Reveive the answers
                for( int process_id = 0; process_id < num_processes; process_id++ ) {
                    int source = ( dipha::mpi_utils::get_rank() + process_id ) % num_processes;
                    mpi_utils::receive_vector( inner_answers_recv_buffer[ source ], source, mpi_utils::MSG_QUERY_SPARSE_INDICES );
                }

                // .. to fill the values of the map
                for( int process_id = 0; process_id < num_processes; process_id++ ) {
                    assert( inner_queries_send_buffer[ process_id ].size() == inner_answers_recv_buffer[ process_id ].size() );
                    for( int64_t query_idx = 0; query_idx < (int64_t)inner_queries_send_buffer[ process_id ].size( ); query_idx++ )
                    {
                        int64_t key = inner_queries_send_buffer[ process_id ][ query_idx ];
                        int64_t value = inner_answers_recv_buffer[ process_id ][ query_idx ];
                        full_to_sparse_indices_in_co_boundary[ key ] = value;
                    }
                }

                // With the map filled, iterate through the queries once again, and answer them
                std::vector< std::vector< int64_t > > answer_buffer;

                std::vector< data_structures::write_once_array_of_arrays< int64_t > > buffer( mpi_utils::get_num_processes() );
                for( int source = 0; source < num_processes; source++ ) {
                    buffer[ source ].init( queries_buffer[ source ].size() );
                    for( int64_t idx = 0; idx < (int64_t)queries_buffer[ source ].size(); idx++ ) {
                        int64_t full_idx = get_locally_full_from_sparse_index( queries_buffer[ source ][ idx ] );
                        if( dual ) {
                            _get_local_coboundary_full_index( full_idx, boundary, true );
                        } else {
                            _get_local_boundary_full_index( full_idx, boundary );
                        }
                        std::vector< int64_t> sparse_idx_boundaries;
                        for( int bd_idx = 0; bd_idx < boundary.size(); bd_idx++ ) {
                            assert( full_to_sparse_indices_in_co_boundary.find( boundary[ bd_idx ] ) != full_to_sparse_indices_in_co_boundary.end() );
                            assert( full_to_sparse_indices_in_co_boundary[ boundary[ bd_idx ] ] != -1 );
                            sparse_idx_boundaries.push_back( full_to_sparse_indices_in_co_boundary[ boundary[ bd_idx ] ] );
                            //std::cout << "add " << boundary[ bd_idx ] << " " << full_to_sparse_indices_in_co_boundary[ boundary[ bd_idx ] ] << std::endl;
                        }
                        buffer[ source ].set( idx, sparse_idx_boundaries.begin(), sparse_idx_boundaries.end() );
                    }
                }

                // send results back to other ranks (including myself)
                std::vector< MPI_Request > answers_requests;
                for( int target = 0; target < num_processes; target++ )
                    mpi_utils::non_blocking_send_vector( buffer[ target ].data, target, mpi_utils::MSG_CO_BOUNDARIES_ANSWERS, answers_requests );

                // receive answers from other ranks
                std::vector< data_structures::write_once_array_of_arrays< int64_t > > answers_of_ranks( num_processes );
                for( int idx = 0; idx < num_processes; idx++ ) {
                    int source = ( mpi_utils::get_rank() + idx ) % num_processes;
                    mpi_utils::receive_vector( answers_of_ranks[ source ].data, source, mpi_utils::MSG_CO_BOUNDARIES_ANSWERS );
                }

                std::vector< int64_t > answers_of_ranks_indices( num_processes, 0 );
                for( int64_t idx = 0; idx < (int64_t)queries.size(); idx++ ) {
                    int target_rank = element_distribution::get_rank( _m_num_elements, queries[ idx ] );
                    const auto& begin = answers_of_ranks[ target_rank ].begin( answers_of_ranks_indices[ target_rank ] );
                    const auto& end = answers_of_ranks[ target_rank ].end( answers_of_ranks_indices[ target_rank ] );
                    answers_of_ranks_indices[ target_rank ]++;
                    answers.set( idx, begin, end );
                }

                // need to make sure that the above sends completed before their buffer goes out of scope
                MPI_Waitall( (int)answers_requests.size(), answers_requests.data(), MPI_STATUSES_IGNORE );

                MPI_Barrier( MPI_COMM_WORLD );

                // Debug: Report result
                /*
                        std::cout << mpi_utils::get_rank() << " Queries ";
                        for( int i = 0; i < queries.size(); i++ ) {
                        std::cout << queries[ i ] << " ";
                        }
                        std::cout << "Answer array: ";
                        for( int j = 0; j < answers.data.size(); j++ ) {
                        std::cout << answers.data[ j ] << " ";
                        }
                        std::cout << std::endl;
                        */

            }

            void _get_local_boundary_full_index( int64_t idx, std::vector< int64_t >& boundary ) const
            {

                boundary.clear();


                std::vector<int64_t> indices;

                int64_t k = _get_local_dim_full_index( idx );

                //std::cout << "k=" << k << std::endl;

                if( k == 0 ) {
                    return;
                }

                idx -= _m_breakpoints[ k ];

                int64_t bcoeff = _m_no_points - 1;
                for( ; k >= 0; k-- ) {
                    while( _m_binomials[ bcoeff ][ k + 1 ] > idx ) {
                        bcoeff--;
                    }

                    indices.push_back( bcoeff );

                    //std::cout << "bcoeff = " << bcoeff << std::endl;

                    idx -= _m_binomials[ bcoeff ][ k + 1 ];

                    bcoeff--;
                }
                std::vector<int64_t>::iterator it = indices.begin(), it2;

                while( it != indices.end() ) {
                    int64_t new_idx;
                    conversion_with_skip( indices.begin(), indices.end(), it, new_idx );
                    boundary.push_back( new_idx );
                    it++;
                }
                return;
            }

            /* Not possible  to determine locally (at least in this way)
               However, we could compute it during load binary and send it around (not done currently)
               void _get_local_coboundary( int64_t idx, std::vector< int64_t >& coboundary ) const
               {
               // Special case, as we have already computed that
               int process_id = dipha::mpi_utils::get_rank();
               assert( _m_breakpoints_local_indices[process_id] <= idx );
               assert( idx < _m_breakpoints_local_indices[process_id+1] );
               idx -= _m_breakpoints_local_indices[process_id];
               coboundary.clear();
               std::copy(_m_coboundaries[idx].begin(), _m_coboundaries[idx].end(), std::back_inserter(coboundary));
               }
               */

        protected:        





            void _precompute_breakpoints()
            {
                _m_breakpoints.push_back( 0 );
                for( int64_t i = 1; i <= _m_upper_dim + 1; i++ ) {
                    _m_breakpoints.push_back( _m_breakpoints[ i - 1 ] + _m_binomials[ _m_no_points ][ i ] );
                }
                /*
                std::cout << "---" << std::endl;
                for(int i=0;i<_m_breakpoints.size();i++) {
                std::cout << _m_breakpoints[i] << std::endl;
                }
                std::cout << "---" << std::endl;
                */

            }

            void _precompute_binomials()
            {
                for( int64_t i = 0; i <= _m_no_points; i++ ) {
                    std::vector<int64_t> binom_row;
                    for( int64_t j = 0; j <= _m_upper_dim + 1; j++ ) {
                        binom_row.push_back( binom( i, j ) );
                    }
                    _m_binomials.push_back( binom_row );
                }
            }

            // Very naive computation, todo: boost (or something else)
            int64_t binom( int64_t n, int64_t k ) const
            {
                if( n < k ) {
                    return 0;
                }
                int64_t result = 1;
                for( int64_t i = n - k + 1; i <= n; i++ ) {
                    result *= i;
                }
                for( int64_t i = k; i >= 1; i-- ) {
                    result /= i;
                }
                return result;
            }


            // other internal functions



            // index -> set of points (in decreasing index order)
            template<typename OutputIterator>
            OutputIterator conversion( int64_t idx, OutputIterator out ) const
            {
                int64_t k = _get_local_dim_full_index( idx );
                idx -= _m_breakpoints[ k ];
                int64_t bcoeff = _m_no_points - 1;
                for( ; k >= 0; k-- ) {
                    while( _m_binomials[ bcoeff ][ k + 1 ] > idx ) {
                        bcoeff--;
                    }
                    *out++ = bcoeff;
                    idx -= _m_binomials[ bcoeff ][ k + 1 ];
                    bcoeff--;
                }
                return out;
            }



            template<typename InputIterator>
            double diameter( InputIterator begin, InputIterator end ) const
            {
                if( begin == end ) {
                    return 0.;
                }
                double max = 0.;
                InputIterator curr = begin;
                do {
                    for( InputIterator run = curr + 1; run != end; run++ ) {
                        double cdist = _m_distance_matrix[ *curr ][ *run ];
                        if( cdist > max ) {
                            max = cdist;
                        }
                    }
                    curr++;
                } while( curr != end );
                return max / 2;
            }

            // Assumes that [begin,end) is sorted in decreasing order!
            template<typename InputIterator>
            int64_t conversion_with_skip( InputIterator begin, InputIterator end, InputIterator skip, int64_t& ind ) const
            {
                int64_t dist = std::distance( begin, end );
                if( skip != end ) {
                    dist--;
                }
                assert( dist <= (int)_m_upper_dim );
                InputIterator it = begin;
                ind = 0;
                for( int64_t k = dist; k >= 1; k-- ) {
                    if( it == skip ) {
                        it++;
                    }
                    ind += _m_binomials[ *it++ ][ k ];
                }
                ind += _m_breakpoints[ dist - 1 ];
                return ind;
            }

            // Assumes that [begin,end) is sorted in decreasing order!
            template<typename InputIterator>
            int64_t conversion_with_extra_index( InputIterator begin, InputIterator end, int64_t extra, int64_t& ind ) const
            {
                int64_t dist = std::distance( begin, end );
                assert( dist <= (int)_m_upper_dim );
                InputIterator it = begin;
                ind = 0;
                bool handled_extra = false;
                for( int64_t k = dist + 1; k >= 1; k-- ) {
                    if( !handled_extra && ( it == end || extra > *it ) ) {
                        ind += _m_binomials[ extra ][ k ];
                        handled_extra = true;
                    } else {
                        ind += _m_binomials[ *it++ ][ k ];
                    }
                }
                ind += _m_breakpoints[ dist ];
                return ind;
            }

        };

    }
}


