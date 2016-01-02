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
    class write_once_column_array
    {

    protected:
      std::vector< int64_t > column_begin;
      std::vector< int64_t > column_end;
      std::vector< int64_t > cell_index;
      std::vector< int64_t > flat_column_entries;

      int64_t local_begin;
      int64_t local_end;

    public:
      void init(int64_t global_num_cols)
      {
        local_begin = element_distribution::get_local_begin(global_num_cols);
        local_end = element_distribution::get_local_end(global_num_cols);
        const int64_t local_num_cols = local_end - local_begin + 1;
        column_begin.resize(local_num_cols, -1);
        column_end.resize(local_num_cols, -1);
        cell_index.resize(local_num_cols, -1);
      }

      bool empty(int64_t idx) const
      {
        assert(idx >= local_begin && idx < local_end);
        return column_begin[idx - local_begin] == column_end[idx - local_begin];
      }

      std::vector< int64_t >::const_iterator begin(int64_t idx) const {
        return flat_column_entries.cbegin() + column_begin[idx - local_begin]; 
      }

      std::vector< int64_t >::const_iterator end(int64_t idx) const { 
        return flat_column_entries.cbegin() + column_end[idx - local_begin]; };

      const int64_t& back(int64_t idx) const { return *(end(idx) - 1); }
      const int64_t& front(int64_t idx) const { return *begin(idx); }

      template< class InputIterator >
      void set(int64_t idx, InputIterator first, InputIterator last, int64_t cell_idx)
      {
        column_begin[idx - local_begin] = flat_column_entries.size();
        std::copy(first, last, std::back_inserter(flat_column_entries));
        column_end[idx - local_begin] = flat_column_entries.size();
        cell_index[idx - local_begin] = cell_idx;
      }

      const int64_t get_cell_index(int64_t idx) const { return cell_index[idx - local_begin]; }
    };
  }
}