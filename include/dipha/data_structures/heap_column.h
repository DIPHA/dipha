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
    class heap_column : public std::vector< int64_t >
    {

    private:
      std::vector< int64_t > temp_col;
      int64_t pushes_since_last_prune;

    public:
      heap_column() : pushes_since_last_prune(0) {};

      template< typename Iterator >
      void add_column(Iterator first, Iterator last)
      {
        for (Iterator it = first; it != last; it++)
          push(*it);
        pushes_since_last_prune += std::distance(first, last);
        if (2 * pushes_since_last_prune > (int64_t)size() + 1024 * 1024)
          prune();
      }

      void prune()
      {
        if (pushes_since_last_prune > 0)
        {
          temp_col.clear();
          int64_t max_index = pop_max_index();
          while (max_index != -1)
          {
            temp_col.push_back(max_index);
            max_index = pop_max_index();
          }
          for (const auto& index : temp_col)
            push(index);
          pushes_since_last_prune = 0;
        }
      }

      void clear()
      {
        pushes_since_last_prune = 0;
        static_cast<std::vector< int64_t >*>(this)->clear();
      }


      int64_t get_max_index()
      {
        int64_t max_index = pop_max_index();
        if (max_index != -1)
          push(max_index);
        return max_index;
      }

    private:
      void push(int64_t index)
      {
        push_back(index);
        std::push_heap(begin(), end());
      }

      int64_t pop()
      {
        int64_t top = front();
        std::pop_heap(begin(), end());
        pop_back();
        return top;
      }

      int64_t pop_max_index()
      {
        if (empty())
          return -1;
        else
        {
          int64_t max_element = pop();
          while (!empty() && front() == max_element)
          {
            pop();
            if (empty())
              return -1;
            else
              max_element = pop();
          }
          return max_element;
        }
      }
    };
  }
}