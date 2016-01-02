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
    template< typename T >
    class write_once_array_of_arrays
    {
      // we make the internal data representation public so that we can easily transmit the object using MPI
    public:
      std::vector< T > data;

    public:
      void init(int64_t num_arrays)
      {
        // data layout: num_arrays % arrays_begin % arrays_end % flat_array_entries
        data.resize(2 * num_arrays + 1, -1);
        data.front() = num_arrays;
      }

      void clear() { data.clear(); }

      int64_t size() const { return data.front(); }

      typename std::vector< T >::const_iterator begin(int64_t idx) const { return data.cbegin() + get_array_begin(idx); }
      typename std::vector< T >::const_iterator end(int64_t idx) const { return data.cbegin() + get_array_end(idx); }

      template< class InputIterator >
      void set(int64_t idx, InputIterator first, InputIterator last)
      {
        data[idx + 1] = data.size();
        data.insert(data.end(), first, last);
        data[size() + idx + 1] = data.size();
      }

    protected:
      int64_t get_array_begin(int64_t idx) const { return data[idx + 1]; }
      int64_t get_array_end(int64_t idx) const { return data[size() + idx + 1]; }
    };
  }
}