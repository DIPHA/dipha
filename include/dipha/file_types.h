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
  namespace file_types
  {
    enum file_type : int64_t
    {
      DIPHA = 8067171840, // == ( int64_t )'D' * ( int64_t )'i' * ( int64_t )'p' * ( int64_t )'h' * ( int64_t )'a'
      WEIGHTED_BOUNDARY_MATRIX = 0,
      IMAGE_DATA = 1,
      PERSISTENCE_DIAGRAM = 2,
      DISTANCE_MATRIX = 7,
      SPARSE_DISTANCE_MATRIX = 8
    };

    inline bool is_dipha_file(const std::string& filename)
    {
      MPI_File file = mpi_utils::file_open_read_only(filename);
      int64_t first_int64_t;
      mpi_utils::file_read_at(file, 0, first_int64_t);
      MPI_File_close(&file);
      return first_int64_t == DIPHA;
    }

    inline void assert_dipha_type(const std::string& filename)
    {
      if (!is_dipha_file(filename))
      {
        mpi_utils::error_printer_if_root() << filename << " is not a proper DIPHA file (first int64_t does not match magic number)" 
                                           << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
    }

    inline file_type get_file_type(const std::string& filename)
    {
      assert_dipha_type(filename);
      MPI_File file = mpi_utils::file_open_read_only(filename);
      int64_t second_int64_t;
      mpi_utils::file_read_at(file, sizeof(int64_t), second_int64_t);
      MPI_File_close(&file);

      return (file_type)second_int64_t;
    }
  }
}
