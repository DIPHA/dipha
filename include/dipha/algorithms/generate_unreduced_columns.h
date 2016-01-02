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
    template< typename Complex >
    void generate_unreduced_columns(const inputs::abstract_weighted_cell_complex< Complex >& complex,
                                    const data_structures::distributed_vector< int64_t >& filtration_to_cell_map,
                                    const data_structures::distributed_vector< int64_t >& cell_to_filtration_map,
                                    int64_t dim,
                                    bool dualize,
                                    data_structures::flat_column_stack& unreduced_columns)
    {
      const int64_t global_num_cells = complex.get_num_cells();
      const int64_t local_begin = element_distribution::get_local_begin(global_num_cells);
      const int64_t local_end = element_distribution::get_local_end(global_num_cells);
      const int64_t local_num_cells = local_end - local_begin + 1;

      const int64_t num_chunks = (local_num_cells + globals::DIPHA_BLOCK_SIZE - 1) / globals::DIPHA_BLOCK_SIZE;

      std::vector< int64_t > dim_queries;
      std::vector< int64_t > dim_answers;
      std::vector< int64_t > cell_to_filtration_map_queries;
      std::vector< int64_t > cell_to_filtration_map_answers;
      data_structures::heap_column boundary;
      std::vector< int64_t > boundary_queries;
      data_structures::write_once_array_of_arrays< int64_t > boundary_answers;
      for (int64_t chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++)
      {
        const int64_t chunk_begin = (local_end - 1) - chunk_idx * globals::DIPHA_BLOCK_SIZE;
        const int64_t chunk_end = std::max((local_end - 1) - (chunk_idx + 1) * globals::DIPHA_BLOCK_SIZE, local_begin);

        dim_queries.clear();
        for (int64_t idx = chunk_begin; idx >= chunk_end; idx--)
        {
          int64_t cur_cell = filtration_to_cell_map.get_local_value(idx);
          dim_queries.push_back(cur_cell);
        }
        complex.get_global_dims(dim_queries, dim_answers);

        auto dim_answers_iterator = dim_answers.cbegin();
        boundary_queries.clear();
        for (int64_t idx = chunk_begin; idx >= chunk_end; idx--)
        {
          int64_t cur_cell = filtration_to_cell_map.get_local_value(idx);
          if (*dim_answers_iterator++ == dim)
            boundary_queries.push_back(cur_cell);
        }
        if (dualize)
          complex.get_global_coboundaries(boundary_queries, boundary_answers);
        else
          complex.get_global_boundaries(boundary_queries, boundary_answers);

        dim_answers_iterator = dim_answers.cbegin();
        int64_t boundary_answers_index = 0;
        cell_to_filtration_map_queries.clear();
        for (int64_t idx = chunk_begin; idx >= chunk_end; idx--)
        {
          if (*dim_answers_iterator++ == dim)
          {
            boundary.assign(boundary_answers.begin(boundary_answers_index), boundary_answers.end(boundary_answers_index));
            boundary_answers_index++;
            for (const auto& boundary_elem : boundary)
              cell_to_filtration_map_queries.push_back(boundary_elem);
          }
        }
        cell_to_filtration_map.get_global_values(cell_to_filtration_map_queries, cell_to_filtration_map_answers);

        boundary_answers_index = 0;
        dim_answers_iterator = dim_answers.cbegin();
        auto cell_to_filtration_map_answers_iterators = cell_to_filtration_map_answers.cbegin();
        for (int64_t idx = chunk_begin; idx >= chunk_end; idx--)
        {
          if (*dim_answers_iterator++ == dim)
          {
            boundary.assign(boundary_answers.begin(boundary_answers_index), boundary_answers.end(boundary_answers_index));
            boundary_answers_index++;
            for (auto& boundary_elem : boundary)
              boundary_elem = *cell_to_filtration_map_answers_iterators++;

            std::make_heap(boundary.begin(), boundary.end());
            unreduced_columns.push(idx, boundary);
          }
        }
      }
    }
  }
}