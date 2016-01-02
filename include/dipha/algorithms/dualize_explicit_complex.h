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
  namespace algorithms
  {
    inline void dualize_explicit_complex(const std::string& input_filename, const std::string& output_filename)
    {
      MPI_Barrier(MPI_COMM_WORLD);

      if (dipha::file_types::get_file_type(input_filename) != dipha::file_types::WEIGHTED_BOUNDARY_MATRIX)
      {
        mpi_utils::error_printer_if_root() << input_filename << " is not a WEIGHTED_BOUNDARY_MATRIX!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }

      std::vector<int64_t> preamble;
      std::vector< int64_t > dims;
      std::vector< double > values;
      std::vector< int64_t > offsets;
      std::vector< int64_t > entries;
      int64_t global_num_entries;
      MPI_Offset file_size;
      {
        MPI_File input_file = mpi_utils::file_open_read_only(input_filename);


        MPI_File_get_size(input_file, &file_size);

        // read preamble
        mpi_utils::file_read_at_vector(input_file, 0, 5, preamble);

        const int64_t global_num_cells = preamble[3];
        const int64_t local_begin = element_distribution::get_local_begin(global_num_cells);
        const int64_t local_end = element_distribution::get_local_end(global_num_cells);
        const int64_t num_local_cells = local_end - local_begin;

        // read dimension data
        MPI_Offset dimensions_begin = (preamble.size() + local_begin) * sizeof(int64_t);
        mpi_utils::file_read_at_vector(input_file, dimensions_begin, num_local_cells, dims);

        // read values data
        MPI_Offset values_begin = (preamble.size() + global_num_cells + local_begin) * sizeof(int64_t);
        mpi_utils::file_read_at_vector(input_file, values_begin, num_local_cells, values);

        // read offsets data
        MPI_Offset offsets_begin = (preamble.size() + 2 * global_num_cells + local_begin) * sizeof(int64_t);
        mpi_utils::file_read_at_vector(input_file, offsets_begin, num_local_cells + 1, offsets);

        // read global_num_entries
        MPI_Offset num_entries_begin = (preamble.size() + 3 * global_num_cells) * sizeof(int64_t);
        mpi_utils::file_read_at(input_file, num_entries_begin, global_num_entries);

        // read entries data
        int64_t entries_start = offsets.front();
        int64_t num_entries = offsets.back() - entries_start;
        MPI_Offset entries_begin = (preamble.size() + 3 * global_num_cells + 1 + entries_start) * sizeof(int64_t);
        mpi_utils::file_read_at_vector(input_file, entries_begin, num_entries, entries);

        MPI_File_close(&input_file);
      }

      std::vector< int64_t > dual_offsets;
      std::vector< int64_t > dual_entries;
      {
        const int64_t global_num_cells = preamble[3];
        const int64_t local_begin = element_distribution::get_local_begin(global_num_cells);
        const int64_t local_end = element_distribution::get_local_end(global_num_cells);
        const int64_t num_local_cells = local_end - local_begin;

        // generate (row, column) pairs
        std::vector< std::pair< int64_t, int64_t > > idx_pairs;
        for (int64_t col = local_begin; col < local_end; col++)
        {
          for (int64_t row_idx = offsets[col - local_begin]; row_idx < offsets[col + 1 - local_begin]; row_idx++)
          {
            idx_pairs.push_back(std::make_pair(entries[row_idx - offsets.front()], col));
          }
        }

        std::vector< std::vector< std::pair< int64_t, int64_t > > > idx_pairs_per_rank(mpi_utils::get_num_processes());
        for (const auto& idx_pair : idx_pairs)
        {
          int target_rank = element_distribution::get_rank(global_num_cells, idx_pair.first);
          idx_pairs_per_rank[target_rank].push_back(idx_pair);
        }

        std::vector< MPI_Request > requests;
        for (int target = 0; target < mpi_utils::get_num_processes(); target++)
          mpi_utils::non_blocking_send_vector(idx_pairs_per_rank[target], target, mpi_utils::MSG_DUALIZED_COLS, requests);


        std::vector< std::vector< std::pair< int64_t, int64_t > > > idx_pairs_per_rank_buffer(mpi_utils::get_num_processes());
        for (int idx = 0; idx < mpi_utils::get_num_processes(); idx++)
        {
          int source = (mpi_utils::get_rank() + idx) % mpi_utils::get_num_processes();
          mpi_utils::receive_vector(idx_pairs_per_rank_buffer[source], source, mpi_utils::MSG_DUALIZED_COLS);
        }

        std::vector< std::pair< int64_t, int64_t > > new_idx_pairs;
        for (int source = 0; source < mpi_utils::get_num_processes(); source++)
        {
          std::copy(idx_pairs_per_rank_buffer[source].begin(), idx_pairs_per_rank_buffer[source].end(), std::back_inserter(new_idx_pairs));
        }

        std::sort(new_idx_pairs.begin(), new_idx_pairs.end());

        int64_t local_num_pairs = new_idx_pairs.size();
        std::vector< int64_t > local_num_pairs_per_rank;
        mpi_utils::all_gather(local_num_pairs, local_num_pairs_per_rank);
        std::vector< int64_t > cum_sum_local_num_pairs(mpi_utils::get_num_processes() + 1);
        cum_sum_local_num_pairs.front() = 0;
        std::partial_sum(local_num_pairs_per_rank.begin(), local_num_pairs_per_rank.end(), cum_sum_local_num_pairs.begin() + 1);

        int64_t local_entries_begin = cum_sum_local_num_pairs[mpi_utils::get_rank()];
        int64_t local_entries_end = cum_sum_local_num_pairs[mpi_utils::get_rank() + 1];

        dual_offsets.resize(num_local_cells, -1);
        for (int64_t idx = 0; idx < (int64_t)new_idx_pairs.size(); idx++)
        {
          if (dual_offsets[new_idx_pairs[idx].first - local_begin] == -1)
            dual_offsets[new_idx_pairs[idx].first - local_begin] = local_entries_begin + dual_entries.size();
          dual_entries.push_back(new_idx_pairs[idx].second);
          if (new_idx_pairs[idx].first - local_begin + 1 < (int64_t)dual_offsets.size())
            dual_offsets[new_idx_pairs[idx].first - local_begin + 1] = local_entries_begin + dual_entries.size();
        }

        dual_offsets.front() = local_entries_begin;
        for (int64_t idx = 0; idx < (int64_t)dual_offsets.size(); idx++)
        {
          if (dual_offsets[idx] == -1)
          {
            dual_offsets[idx] = dual_offsets[idx - 1];
          }
        }

        MPI_Waitall((int)requests.size(), requests.data(), MPI_STATUSES_IGNORE);
      }

      MPI_File output_file = mpi_utils::file_open_write_only(output_filename);

      MPI_File_set_size(output_file, file_size);

      const auto input_boundary_type = (inputs::weighted_explicit_complex::boundary_type)preamble[2];

      inputs::weighted_explicit_complex::boundary_type output_boundary_type;
      if (input_boundary_type == inputs::weighted_explicit_complex::boundary_type::BOUNDARY)
        output_boundary_type = inputs::weighted_explicit_complex::boundary_type::COBOUNDARY;
      else
        output_boundary_type = inputs::weighted_explicit_complex::boundary_type::BOUNDARY;

      preamble[2] = (int64_t)output_boundary_type;

      if (mpi_utils::is_root())
        mpi_utils::file_write_at_vector(output_file, 0, preamble);

      const int64_t global_num_cells = preamble[3];
      const int64_t local_begin = element_distribution::get_local_begin(global_num_cells);
      const int64_t local_end = element_distribution::get_local_end(global_num_cells);

      MPI_Offset dimensions_begin = (preamble.size() + local_begin) * sizeof(int64_t);
      mpi_utils::file_write_at_vector(output_file, dimensions_begin, dims);

      MPI_Offset values_begin = (preamble.size() + global_num_cells + local_begin) * sizeof(int64_t);
      mpi_utils::file_write_at_vector(output_file, values_begin, values);

      MPI_Offset offsets_begin = (preamble.size() + 2 * global_num_cells + local_begin) * sizeof(int64_t);
      mpi_utils::file_write_at_vector(output_file, offsets_begin, dual_offsets);

      MPI_Offset num_entries_begin = (preamble.size() + 3 * global_num_cells) * sizeof(int64_t);
      if (mpi_utils::is_root())
        mpi_utils::file_write_at(output_file, num_entries_begin, global_num_entries);

      int64_t entries_start = dual_offsets.front();
      MPI_Offset entries_begin = (preamble.size() + 3 * global_num_cells + 1 + entries_start) * sizeof(int64_t);
      mpi_utils::file_write_at_vector(output_file, entries_begin, dual_entries);

      MPI_File_close(&output_file);

      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}