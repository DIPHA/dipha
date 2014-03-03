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
        // see http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern for template voodoo below
        template< typename Derived >
        class abstract_weighted_cell_complex {

        protected:
            const Derived& derived() const { return static_cast<const Derived&>( *this ); }

            Derived& derived() { return static_cast<Derived&>( *this ); }

            // functions that MUST be implemented by template parameter
        public:
            // total number of cells in the complex
            int64_t get_num_cells() const { return derived()._get_num_cells(); }

            // dimension of complex
            int64_t get_max_dim() const { return derived()._get_max_dim(); }

            // dimension of given cell
            int64_t get_local_dim( int64_t idx ) const
            {
                assert( idx >= element_distribution::get_local_begin( get_num_cells() ) && idx < element_distribution::get_local_end( get_num_cells() ) );
                return derived()._get_local_dim( idx );
            }

            // value of given cell
            double get_local_value( int64_t idx ) const
            {
                assert( idx >= element_distribution::get_local_begin( get_num_cells() ) && idx < element_distribution::get_local_end( get_num_cells() ) );
                return derived()._get_local_value( idx );
            }

            // boundary of given cell
            void get_local_boundary( int64_t idx, std::vector< int64_t >& boundary ) const
            {
                assert( idx >= element_distribution::get_local_begin( get_num_cells() ) && idx < element_distribution::get_local_end( get_num_cells() ) );
                boundary.clear();
                derived()._get_local_boundary( idx, boundary );
            }

            // coboundary of given cell
            void get_local_coboundary( int64_t idx, std::vector< int64_t >& coboundary ) const
            {
                assert( idx >= element_distribution::get_local_begin( get_num_cells() ) && idx < element_distribution::get_local_end( get_num_cells() ) );
                coboundary.clear();
                derived()._get_local_coboundary( idx, coboundary );
            }

            // loads the object from file
            void load_binary( const std::string& filename )
            {
                file_types::assert_dipha_type( filename );
                MPI_File file = mpi_utils::file_open_read_only( filename );
                derived()._load_binary( file );
                MPI_File_close( &file );
            }

            // functions that CAN be implemented by template parameter to improve performance
        public:
            void get_global_values( const std::vector< int64_t >& queries,
                                    std::vector< double >& answers ) const
            {
                answers.clear();
                answers.reserve( queries.size() );
                derived()._get_global_values( queries, answers );
            }

            void get_global_dims( const std::vector< int64_t >& queries,
                                  std::vector< int64_t >& answers ) const
            {
                answers.clear();
                answers.reserve( queries.size() );
                derived()._get_global_dims( queries, answers );
            }

            void get_global_boundaries( const std::vector< int64_t >& queries,
                                        data_structures::write_once_array_of_arrays< int64_t >& answers ) const
            {
                answers.clear();
                answers.init( queries.size() );
                derived()._get_global_boundaries( queries, answers );
            }

            void get_global_coboundaries( const std::vector< int64_t >& queries,
                                          data_structures::write_once_array_of_arrays< int64_t >& answers ) const
            {
                answers.clear();
                answers.init( queries.size() );
                derived()._get_global_coboundaries( queries, answers );
            }

            double get_max_value() const
            {
                return derived()._get_max_value();
            }

            // default implementations
        private:
            void _get_global_values( const std::vector< int64_t >& queries,
                                     std::vector< double >& answers ) const
            {

                const int64_t global_num_cells = get_num_cells();
                std::vector< std::vector< int64_t > > queries_buffer;
                element_distribution::scatter_queries( queries, global_num_cells, queries_buffer );

                // process queries
                std::vector< std::vector< double > > answers_buffer( mpi_utils::get_num_processes() );
                for( int source = 0; source < mpi_utils::get_num_processes(); source++ ) {
                    for( const auto& query : queries_buffer[ source ] )
                        answers_buffer[ source ].push_back( get_local_value( query ) );
                }

                element_distribution::gather_answers( queries, global_num_cells, answers_buffer, answers );
            }

            void _get_global_dims( const std::vector< int64_t >& queries,
                                   std::vector< int64_t >& answers ) const
            {
                const int64_t global_num_cells = get_num_cells();
                std::vector< std::vector< int64_t > > queries_buffer;
                element_distribution::scatter_queries( queries, global_num_cells, queries_buffer );

                // process queries
                std::vector< std::vector< int64_t > > answers_buffer( mpi_utils::get_num_processes() );
                for( int source = 0; source < mpi_utils::get_num_processes(); source++ ) {
                    for( const auto& query : queries_buffer[ source ] )
                        answers_buffer[ source ].push_back( get_local_dim( query ) );
                }

                element_distribution::gather_answers( queries, global_num_cells, answers_buffer, answers );
            }

            void _get_global_boundaries( const std::vector< int64_t >& queries,
                                         data_structures::write_once_array_of_arrays< int64_t >& answers ) const
            {
                _get_global_co_boundaries( queries, answers, false );
            }

            void _get_global_coboundaries( const std::vector< int64_t >& queries,
                                           data_structures::write_once_array_of_arrays< int64_t >& answers ) const
            {
                _get_global_co_boundaries( queries, answers, true );
            }

            double _get_max_value() const
            {
                const int64_t local_begin = element_distribution::get_local_begin( get_num_cells( ) );
                const int64_t local_end = element_distribution::get_local_end( get_num_cells( ) );
                double local_max_value = std::numeric_limits< double >::lowest( );
                for( int64_t idx = local_begin; idx < local_end; idx++ ) {
                    double value = get_local_value( idx );
                    local_max_value = value > local_max_value ? value : local_max_value;
                }
                std::vector< double > max_value_per_rank( mpi_utils::get_num_processes( ) );
                MPI_Allgather( &local_max_value, 1, MPI_DOUBLE, max_value_per_rank.data( ), 1, MPI_DOUBLE, MPI_COMM_WORLD );
                return *std::max_element( max_value_per_rank.begin( ), max_value_per_rank.end( ) );
            }

            // internal helper functions
        private:
            void _get_global_co_boundaries( const std::vector< int64_t >& queries,
                                            data_structures::write_once_array_of_arrays< int64_t >& answers,
                                            bool dual ) const
            {
                const int64_t global_num_cells = get_num_cells( );

                std::vector< std::vector< int64_t > > queries_buffer;
                element_distribution::scatter_queries( queries, global_num_cells, queries_buffer );

                // process queries
                std::vector< int64_t > boundary;
                std::vector< data_structures::write_once_array_of_arrays< int64_t > > buffer( mpi_utils::get_num_processes( ) );
                for( int source = 0; source < mpi_utils::get_num_processes( ); source++ ) {
                    buffer[ source ].init( queries_buffer[ source ].size( ) );
                    for( int64_t idx = 0; idx < (int64_t)queries_buffer[ source ].size( ); idx++ ) {
                        if( dual )
                            get_local_coboundary( queries_buffer[ source ][ idx ], boundary );
                        else
                            get_local_boundary( queries_buffer[ source ][ idx ], boundary );
                        buffer[ source ].set( idx, boundary.begin( ), boundary.end( ) );
                    }
                }

                // send results back to other ranks (including myself)
                std::vector< MPI_Request > answers_requests;
                for( int target = 0; target < mpi_utils::get_num_processes( ); target++ )
                    mpi_utils::non_blocking_send_vector( buffer[ target ].data, target, mpi_utils::MSG_CO_BOUNDARIES_ANSWERS, answers_requests );

                // receive answers from other ranks
                std::vector< data_structures::write_once_array_of_arrays< int64_t > > answers_of_ranks( mpi_utils::get_num_processes( ) );
                for( int idx = 0; idx < mpi_utils::get_num_processes( ); idx++ ) {
                    int source = ( mpi_utils::get_rank( ) + idx ) % mpi_utils::get_num_processes( );
                    mpi_utils::receive_vector( answers_of_ranks[ source ].data, source, mpi_utils::MSG_CO_BOUNDARIES_ANSWERS );
                }

                std::vector< int64_t > answers_of_ranks_indices( mpi_utils::get_num_processes( ), 0 );
                for( int64_t idx = 0; idx < (int64_t)queries.size( ); idx++ ) {
                    int target_rank = element_distribution::get_rank( global_num_cells, queries[ idx ] );
                    const auto& begin = answers_of_ranks[ target_rank ].begin( answers_of_ranks_indices[ target_rank ] );
                    const auto& end = answers_of_ranks[ target_rank ].end( answers_of_ranks_indices[ target_rank ] );
                    answers_of_ranks_indices[ target_rank ]++;
                    answers.set( idx, begin, end );
                }

                // need to make sure that the above sends completed before their buffer goes out of scope
                MPI_Waitall( (int)answers_requests.size( ), answers_requests.data( ), MPI_STATUSES_IGNORE );

                MPI_Barrier( MPI_COMM_WORLD );
            }
        };
    }
}

#ifdef DIPHA_TEST
    template< typename T >
    class TestPrimalWeightedCellComplex : public ::testing::Test {
        virtual void SetUp() {
            this->complex.load_binary( T::get_test_filename() );
        }
    public:
        typename T::complex_type complex;
    };
    TYPED_TEST_CASE_P( TestPrimalWeightedCellComplex );

    TYPED_TEST_P( TestPrimalWeightedCellComplex, BoundaryDim )
    {
        const int64_t num_cells = this->complex.get_num_cells();

        std::vector< int64_t > cells;
        dipha::data_structures::write_once_array_of_arrays< int64_t > boundaries;
        std::vector< int64_t > cells_dim;
        for( int64_t idx = 0; idx < num_cells; idx++ )
            cells.push_back( idx );
        this->complex.get_global_boundaries( cells, boundaries );
        this->complex.get_global_dims( cells, cells_dim );

        std::vector< int64_t > boundaries_entries;
        std::vector< int64_t > boundaries_entries_dim;
        for( int64_t idx = 0; idx < num_cells; idx++ ) {
            for( auto it = boundaries.begin( idx ); it != boundaries.end( idx ); it++ ) {
                boundaries_entries.push_back( *it );
            }
        }
        this->complex.get_global_dims( boundaries_entries, boundaries_entries_dim );

        auto iterator_of_boundaries_entries_dim = boundaries_entries_dim.cbegin( );
        auto iterator_of_cells_dim = cells_dim.cbegin( );

        for( int64_t idx = 0; idx < num_cells; idx++ ) {
            int64_t dim_of_cell = *iterator_of_cells_dim++;
            for( auto it = boundaries.begin( idx ); it != boundaries.end( idx ); it++ ) {
                int64_t dim_of_boundary_entry = *iterator_of_boundaries_entries_dim++;
                ASSERT_EQ( dim_of_boundary_entry, dim_of_cell - 1 );
            }
        }
    }

    TYPED_TEST_P( TestPrimalWeightedCellComplex, AnotherPrimalTest )
    {
        //...
    }

    REGISTER_TYPED_TEST_CASE_P( TestPrimalWeightedCellComplex, BoundaryDim, AnotherPrimalTest );


    template< typename T >
    class TestDualWeightedCellComplex : public ::testing::Test {
        virtual void SetUp( )
        {
            this->complex.load_binary( T::get_test_filename( ) );
        }
    public:
        typename T::complex_type complex;
    };
    TYPED_TEST_CASE_P( TestDualWeightedCellComplex );

    TYPED_TEST_P( TestDualWeightedCellComplex, CoboundaryDim )
    {
        const int64_t num_cells = this->complex.get_num_cells( );

        std::vector< int64_t > cells;
        dipha::data_structures::write_once_array_of_arrays< int64_t > coboundaries;
        std::vector< int64_t > cells_dim;
        for( int64_t idx = 0; idx < num_cells; idx++ )
            cells.push_back( idx );
        this->complex.get_global_coboundaries( cells, coboundaries );
        this->complex.get_global_dims( cells, cells_dim );

        std::vector< int64_t > coboundaries_entries;
        std::vector< int64_t > coboundaries_entries_dim;
        for( int64_t idx = 0; idx < num_cells; idx++ ) {
            for( auto it = coboundaries.begin( idx ); it != coboundaries.end( idx ); it++ ) {
                coboundaries_entries.push_back( *it );
            }
        }
        this->complex.get_global_dims( coboundaries_entries, coboundaries_entries_dim );

        auto iterator_of_coboundaries_entries_dim = coboundaries_entries_dim.cbegin( );
        auto iterator_of_cells_dim = cells_dim.cbegin( );

        for( int64_t idx = 0; idx < num_cells; idx++ ) {
            int64_t dim_of_cell = *iterator_of_cells_dim++;
            for( auto it = coboundaries.begin( idx ); it != coboundaries.end( idx ); it++ ) {
                int64_t dim_of_coboundary_entry = *iterator_of_coboundaries_entries_dim++;
                ASSERT_EQ( dim_of_coboundary_entry, dim_of_cell + 1 );
            }
        }
    }

    TYPED_TEST_P( TestDualWeightedCellComplex, AnotherDualTest )
    {
        //...
    }

    REGISTER_TYPED_TEST_CASE_P( TestDualWeightedCellComplex, CoboundaryDim, AnotherDualTest );
#endif


