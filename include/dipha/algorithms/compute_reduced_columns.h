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
    template<typename Complex>
    void compute_reduced_columns(const inputs::abstract_weighted_cell_complex<Complex>& complex,
                                 bool dualize,
                                 int64_t upper_dim,
                                 data_structures::distributed_vector<int64_t>& filtration_to_cell_map,
                                 data_structures::write_once_column_array& reduced_columns)
    {
      DIPHA_MACROS_BENCHMARK(get_filtration_to_cell_map(complex, dualize, filtration_to_cell_map); );
      data_structures::distributed_vector<int64_t> cell_to_filtration_map;
      DIPHA_MACROS_BENCHMARK(get_cell_to_filtration_map(complex.get_num_cells(), filtration_to_cell_map, cell_to_filtration_map); );

      reduced_columns.init(complex.get_num_cells());
      const int64_t max_dim = std::min(upper_dim, complex.get_max_dim());
      for (int64_t idx = 0; idx < max_dim; idx++)
      {
        int64_t cur_dim = dualize ? idx : max_dim - idx;
        data_structures::flat_column_stack unreduced_columns;
        DIPHA_MACROS_BENCHMARK(generate_unreduced_columns(complex, filtration_to_cell_map, cell_to_filtration_map, 
                                                          cur_dim, dualize, unreduced_columns); );
        DIPHA_MACROS_BENCHMARK(reduction_kernel(complex.get_num_cells(), unreduced_columns, reduced_columns); );
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
}