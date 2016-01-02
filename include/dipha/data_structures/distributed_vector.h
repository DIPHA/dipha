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
  namespace data_structures
  {
    template< typename ValueType >
    class distributed_vector
    {

    public:
      typedef ValueType value_type;

    protected:
      std::vector< value_type > entries;
      int64_t local_begin;
      int64_t local_end;
      int64_t global_num_elems;

    public:
      void init(int64_t arg_global_num_elems)
      {
        local_begin = element_distribution::get_local_begin(arg_global_num_elems);
        local_end = element_distribution::get_local_end(arg_global_num_elems);
        const int64_t local_num_elems = local_end - local_begin;
        entries.resize(local_num_elems);
        global_num_elems = arg_global_num_elems;
      }

      void set_local_value(int64_t idx, const value_type& value)
      {
        assert(idx >= local_begin && idx < local_end);
        entries[idx - local_begin] = value;
      }

      value_type get_local_value(int64_t idx) const
      {
        assert(idx >= local_begin && idx < local_end);
        return entries[idx - local_begin];
      }

      void get_global_values(const std::vector< int64_t >& indices,
                             std::vector< value_type >& values) const
      {
        std::vector< std::vector< int64_t > > queries_buffer;
        element_distribution::scatter_queries(indices, global_num_elems, queries_buffer);

        // process queries
        std::vector< std::vector< int64_t > > answers_buffer(mpi_utils::get_num_processes());
        for (int source = 0; source < mpi_utils::get_num_processes(); source++)
        {
          for (const auto& query : queries_buffer[source])
            answers_buffer[source].push_back(entries[query - local_begin]);
        }

        element_distribution::gather_answers(indices, global_num_elems, answers_buffer, values);
      }


      void set_global_values(const std::vector< std::pair< int64_t, value_type > >& indices_and_values)
      {
        // 1. generate pairs for all ranks
        std::vector< std::vector< std::pair< int64_t, value_type > > > pairs_for_ranks(mpi_utils::get_num_processes());
        for (const auto& index_and_value : indices_and_values)
          pairs_for_ranks[element_distribution::get_rank(global_num_elems, index_and_value.first)].push_back(index_and_value);

        // 2. send pairs to other ranks
        std::vector< MPI_Request > requests;
        for (int target = 0; target < mpi_utils::get_num_processes(); target++)
          mpi_utils::non_blocking_send_vector(pairs_for_ranks[target], target, mpi_utils::MSG_SET_GLOBAL_VALUES, requests);

        // 3. receive pairs from other ranks
        std::vector< std::pair< int64_t, value_type > > buffer;
        for (int idx = 0; idx < mpi_utils::get_num_processes(); idx++)
        {
          const int source = (mpi_utils::get_rank() + idx) % mpi_utils::get_num_processes();
          mpi_utils::receive_vector(buffer, source, mpi_utils::MSG_SET_GLOBAL_VALUES);
          for (const auto& buffer_elem : buffer)
            set_local_value(buffer_elem.first, buffer_elem.second);
        }

        // 4. need to make sure that above send completed before its buffer goes out of scope
        MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
      }
    };
  }
}