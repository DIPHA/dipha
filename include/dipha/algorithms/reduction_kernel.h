/*  Copyright 2014 IST Austria

Contributed by: Jan Reininghaus, Michael Kerber, Ulrich Bauer

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
  namespace algorithms
  {
    inline void reduction_kernel(int64_t global_num_cols,
                                 data_structures::flat_column_stack& unreduced_columns,
                                 data_structures::write_once_column_array& reduced_columns)
    {
      const int64_t local_begin = element_distribution::get_local_begin(global_num_cols);
      const double reduction_kernel_start = MPI_Wtime();
      const int num_processes = mpi_utils::get_num_processes();
      const int rank = mpi_utils::get_rank();

      data_structures::heap_column col;
      for (int round = 0; round < num_processes - rank; round++)
      {
        data_structures::flat_column_stack temp_columns;
        while (!unreduced_columns.empty())
        {
          int64_t index;
          unreduced_columns.pop(index, col);

          // apply twist optimization
          if (round == 0 && !reduced_columns.empty(index))
            col.clear();

          int64_t maximal_index = col.get_max_index();

          // (partially) reduce col using reduced_columns
          if (maximal_index >= local_begin && !reduced_columns.empty(maximal_index))
          {
            while (maximal_index >= local_begin && !reduced_columns.empty(maximal_index))
            {
              col.add_column(reduced_columns.begin(maximal_index), reduced_columns.end(maximal_index));
              maximal_index = col.get_max_index();
            }
            col.prune();
          }

          if (!col.empty())
          {
            // store locally or pass col on?
            if (maximal_index >= local_begin)
              reduced_columns.set(maximal_index, col.begin(), col.end(), index);
            else
              temp_columns.push(index, col);
          }

          unreduced_columns.shrink_to_fit();
        }

        // reverse column stack
        while (!temp_columns.empty())
        {
          int64_t index;
          temp_columns.pop(index, col);
          unreduced_columns.push(index, col);
        }

        // depending on rank: send unreduced_columns to next rank
        std::vector< MPI_Request > requests;
        if (rank > 0)
          mpi_utils::non_blocking_send_vector(unreduced_columns.stack_data, rank - 1, mpi_utils::MSG_UNREDUCED_COLUMNS, requests);

        // depending on round and rank: receive unreduced_columns from previous rank
        if (round < num_processes - rank - 1)
          mpi_utils::receive_vector(temp_columns.stack_data, rank + 1, mpi_utils::MSG_UNREDUCED_COLUMNS);

        if (rank > 0)
          MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);

        temp_columns.swap(unreduced_columns);
      }

      globals::reduction_kernel_running_time += MPI_Wtime() - reduction_kernel_start;
    }
  }
}