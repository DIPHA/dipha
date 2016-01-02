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

namespace dipha
{
  namespace element_distribution
  {
    inline int64_t get_range_width(int64_t global_num_elems)
    {
      return (global_num_elems + mpi_utils::get_num_processes() - 1) / mpi_utils::get_num_processes();
    }

    inline int get_rank(int64_t global_num_elems, int64_t elem)
    {
      return (int)(elem / get_range_width(global_num_elems));
    }

    inline int64_t get_local_begin(int64_t global_num_elems, int rank)
    {
      return std::min(get_range_width(global_num_elems) * rank, global_num_elems);
    }

    inline int64_t get_local_begin(int64_t global_num_elems)
    {
      return get_local_begin(global_num_elems, mpi_utils::get_rank());
    }

    inline int64_t get_local_end(int64_t global_num_elems, int rank)
    {
      return std::min(get_range_width(global_num_elems) * (rank + 1), global_num_elems);
    }

    inline int64_t get_local_end(int64_t global_num_elems)
    {
      return get_local_end(global_num_elems, mpi_utils::get_rank());
    }

    inline void scatter_queries(const std::vector< int64_t >& queries,
                                int64_t global_num_cells,
                                std::vector< std::vector< int64_t > >& queries_buffer)
    {
      MPI_Barrier(MPI_COMM_WORLD);

      std::vector< std::vector< int64_t > > queries_for_ranks(mpi_utils::get_num_processes());
      for (const auto& query : queries)
        queries_for_ranks[get_rank(global_num_cells, query)].push_back(query);

      // send queries to other ranks (including myself)
      std::vector< MPI_Request > queries_requests;
      for (int target = 0; target < mpi_utils::get_num_processes(); target++)
        mpi_utils::non_blocking_send_vector(queries_for_ranks[target], target, mpi_utils::MSG_SCATTER_QUERIES, queries_requests);

      // receive requests from other ranks
      queries_buffer.resize(mpi_utils::get_num_processes());
      for (int idx = 0; idx < mpi_utils::get_num_processes(); idx++)
      {
        int source = (mpi_utils::get_rank() + idx) % mpi_utils::get_num_processes();
        mpi_utils::receive_vector(queries_buffer[source], source, mpi_utils::MSG_SCATTER_QUERIES);
      }

      MPI_Waitall((int)queries_requests.size(), queries_requests.data(), MPI_STATUSES_IGNORE);
    }

    template< typename T >
    void gather_answers(const std::vector< int64_t >& queries,
                        int64_t global_num_cells,
                        const std::vector< std::vector< T > >& answers_buffer,
                        std::vector< T >& answers)
    {
      MPI_Barrier(MPI_COMM_WORLD);

      // send results back to other ranks (including myself)
      std::vector< MPI_Request > answers_requests;
      for (int target = 0; target < mpi_utils::get_num_processes(); target++)
        mpi_utils::non_blocking_send_vector(answers_buffer[target], target, mpi_utils::MSG_GATHER_ANSWERS, answers_requests);

      // receive answers from other ranks
      std::vector< std::vector< T > > answers_of_ranks(mpi_utils::get_num_processes());
      for (int idx = 0; idx < mpi_utils::get_num_processes(); idx++)
      {
        int source = (mpi_utils::get_rank() + idx) % mpi_utils::get_num_processes();
        mpi_utils::receive_vector(answers_of_ranks[source], source, mpi_utils::MSG_GATHER_ANSWERS);
      }

      std::vector< typename std::vector< T >::const_iterator > answers_of_ranks_iterators;
      for (const auto& answers_of_rank : answers_of_ranks)
        answers_of_ranks_iterators.push_back(answers_of_rank.cbegin());

      answers.clear();
      for (const auto& query : queries)
        answers.push_back(*answers_of_ranks_iterators[get_rank(global_num_cells, query)]++);

      MPI_Waitall((int)answers_requests.size(), answers_requests.data(), MPI_STATUSES_IGNORE);
    }
  }
}
