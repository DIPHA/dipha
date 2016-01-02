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
  namespace inputs
  {
    class weighted_explicit_complex : public abstract_weighted_cell_complex< weighted_explicit_complex >
    {

      friend abstract_weighted_cell_complex< weighted_explicit_complex >;

    public:
      enum class boundary_type : int64_t { BOUNDARY = 0, COBOUNDARY = 1 };

      // defining state of the object
    protected:
      std::vector< int64_t > offsets;
      std::vector< int64_t > entries;
      std::vector< double > values;
      std::vector< int64_t > dims;
      boundary_type my_boundary_type;
      int64_t global_num_cells;
      int64_t max_dim;

      // derived quantities
    protected:
      int64_t local_begin;

      // implementation of abstract_weighted_cell_complex interface
    protected:
      int64_t _get_num_cells() const { return global_num_cells; }

      int64_t _get_local_dim(int64_t idx) const { return dims[idx - local_begin]; }

      double _get_local_value(int64_t idx) const { return values[idx - local_begin]; }

      void _get_local_boundary(int64_t idx, std::vector< int64_t >& boundary) const { _get_local_co_boundary(idx, false, boundary); }

      void _get_local_coboundary(int64_t idx, std::vector< int64_t >& coboundary) const { _get_local_co_boundary(idx, true, coboundary); }

      int64_t _get_max_dim() const { return max_dim; }

      // Loads the weighted_cubical_complex from given file in binary format -- all symbols are 64 bit wide
      // Format: file_types::DIPHA % file_types::WEIGHTED_BOUNDARY_MATRIX % boundary_type % num_cells (N) % max_dim % dim1 % ... 
      //         % dimN % value1 % ... % valueN % offset1 % ... % offsetN % num_entries (M) % entry1 % ... % entryM 
      void _load_binary(MPI_File file,
                        int64_t upper_dim = std::numeric_limits< int64_t >::max())
      {
        // read preamble
        std::vector< int64_t > preamble;
        mpi_utils::file_read_at_vector(file, 0, 5, preamble);
        int64_t dipha_identifier = preamble[0];
        int64_t file_type = preamble[1];
        my_boundary_type = (boundary_type)preamble[2];
        global_num_cells = preamble[3];
        max_dim = preamble[4];

        local_begin = element_distribution::get_local_begin(global_num_cells);
        const int64_t num_local_cells = element_distribution::get_local_end(global_num_cells) - local_begin;

        // read dimension data
        MPI_Offset dimensions_begin = (preamble.size() + local_begin) * sizeof(int64_t);
        mpi_utils::file_read_at_vector(file, dimensions_begin, num_local_cells, dims);

        // read values data
        MPI_Offset values_begin = (preamble.size() + global_num_cells + local_begin) * sizeof(int64_t);
        mpi_utils::file_read_at_vector(file, values_begin, num_local_cells, values);

        // read offsets data
        MPI_Offset offsets_begin = (preamble.size() + 2 * global_num_cells + local_begin) * sizeof(int64_t);
        mpi_utils::file_read_at_vector(file, offsets_begin, num_local_cells + 1, offsets);

        // read entries data
        int64_t entries_start = offsets.front();
        int64_t local_num_entries = offsets.back() - entries_start;
        MPI_Offset entries_begin = (preamble.size() + 3 * global_num_cells + 1 + entries_start) * sizeof(int64_t);
        mpi_utils::file_read_at_vector(file, entries_begin, local_num_entries, entries);
      }

      // internal helper functions
    protected:
      void _get_local_co_boundary(int64_t idx, bool dual, std::vector< int64_t >& co_boundary) const
      {
        if (!dual && my_boundary_type == boundary_type::COBOUNDARY)
        {
          mpi_utils::error_printer_if_root() << "Primal computation not supported for this input file. " 
                                             << "Please convert it first using the dualize utility.";
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        if (dual && my_boundary_type == boundary_type::BOUNDARY)
        {
          mpi_utils::error_printer_if_root() << "Dual computation not supported for this input file. " 
                                             << "Please convert it first using the dualize utility.";
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        co_boundary.clear();
        for (int64_t entry_idx = offsets[idx - local_begin]; entry_idx < offsets[idx + 1 - local_begin]; entry_idx++)
        {
          co_boundary.push_back(entries[entry_idx - offsets.front()]);
        }
      }
    };
  }
}
