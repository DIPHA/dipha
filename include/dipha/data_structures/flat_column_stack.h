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
    class flat_column_stack
    {

      // we make the internal data representation public so that we can easily transmit the object using MPI
    public:
      std::vector< int64_t > stack_data;

    public:
      void swap(flat_column_stack& other) { std::swap(other.stack_data, stack_data); }
      void clear() { stack_data.clear(); }
      bool empty() const { return stack_data.empty(); }

      void shrink_to_fit()
      {
        if (stack_data.capacity() > 2 * stack_data.size())
          stack_data.shrink_to_fit();
      }

      void push(int64_t index, const heap_column& col)
      {
        stack_data.push_back(index);
        stack_data.insert(stack_data.end(), col.begin(), col.end());
        stack_data.push_back(col.size());
      }

      void pop(int64_t& index, heap_column& col)
      {
        const int64_t col_size = stack_data.back();
        col.clear();
        col.resize(col_size);
        std::copy(stack_data.rbegin() + 1, stack_data.rbegin() + 1 + col_size, col.rbegin());
        index = stack_data[stack_data.size() - 2 - col_size];
        stack_data.resize(stack_data.size() - 2 - col_size);
      }
    };
  }
}