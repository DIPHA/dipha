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
    inline void get_cell_to_filtration_map(int64_t global_num_cells,
                                           const data_structures::distributed_vector< int64_t >& filtration_to_cell_map,
                                           data_structures::distributed_vector< int64_t >& cell_to_filtration_map)
    {
      const int64_t local_begin = element_distribution::get_local_begin(global_num_cells);
      const int64_t local_end = element_distribution::get_local_end(global_num_cells);
      const int64_t local_num_cells = local_end - local_begin;

      std::vector < std::pair< int64_t, int64_t > > indices_and_values;
      for (int64_t idx = local_begin; idx < local_end; idx++)
        indices_and_values.push_back(std::make_pair(filtration_to_cell_map.get_local_value(idx), idx));

      cell_to_filtration_map.init(global_num_cells);
      cell_to_filtration_map.set_global_values(indices_and_values);
    }
  }
}