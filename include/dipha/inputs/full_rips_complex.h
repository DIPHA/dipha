/*  Copyright 2014 IST Austria

Contributed by: Michael Kerber

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
    class full_rips_complex : public abstract_weighted_cell_complex< full_rips_complex >
    {

      friend abstract_weighted_cell_complex< full_rips_complex >;

      // defining state of the object
    protected:
      int64_t _m_upper_dim;
      int64_t _m_no_points;

      std::vector<std::vector<double> > _m_distance_matrix;

      // derived quantities
    protected:

      std::vector<std::vector<int64_t> > _m_binomials;

      // contains the indices where a new dimension starts
      std::vector<int64_t> _m_breakpoints;

    protected:

      // Loads the complete_rips_complex from given file in binary format -- all symbols are 64 bit wide
      void _load_binary(MPI_File file,
                        int64_t upper_dim = std::numeric_limits< int64_t >::max())
      {
        // read preamble
        std::vector< int64_t > preamble;
        mpi_utils::file_read_at_vector(file, 0, 3, preamble);

        _m_no_points = preamble[2];

        if (upper_dim == std::numeric_limits< int64_t >::max())
        {
          mpi_utils::error_printer_if_root() << "No upper_dim specified for distance_matrix data!";
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        else
        {
          _m_upper_dim = upper_dim;
        }

        int64_t matrix_size = _m_no_points * _m_no_points;

        std::vector< double > distances;
        MPI_Offset offset = preamble.size() * sizeof(int64_t);
        mpi_utils::file_read_at_vector(file, offset, matrix_size, distances);


        _m_distance_matrix.resize(_m_no_points);
        for (int64_t i = 0; i < _m_no_points; i++)
        {
          _m_distance_matrix[i].resize(_m_no_points);
        }

        for (int64_t i = 0; i < matrix_size; i++)
        {
          _m_distance_matrix[i / _m_no_points][i % _m_no_points] = distances[i];
        }

        _precompute_binomials();

        _precompute_breakpoints();
      }

    public:

      double get_distance(int64_t i, int64_t j) const
      {
        return _m_distance_matrix[i][j];
      }

      int64_t number_of_points() const
      {
        return _m_no_points;
      }

    protected:

      int64_t _get_max_dim() const
      {
        return _m_upper_dim;
      }

      int64_t _get_num_cells() const
      {
        return _m_breakpoints[_m_upper_dim + 1];
      }

      int64_t _get_local_dim(int64_t idx) const
      {
        int64_t k = 0;
        while (k < _m_upper_dim + 1 && _m_breakpoints[k + 1] <= idx)
        {
          k++;
        }
        return k;
      }

      double _get_local_value(int64_t idx) const
      {
        static std::vector<int64_t> points;
        points.clear();
        conversion(idx, std::back_inserter(points));
        return diameter(points.begin(), points.end());
      }


      void _get_local_boundary(int64_t idx, std::vector< int64_t >& boundary) const
      {

        boundary.clear();


        std::vector<int64_t> indices;

        int64_t k = _get_local_dim(idx);

        //std::cout << "k=" << k << std::endl;

        if (k == 0)
        {
          return;
        }

        idx -= _m_breakpoints[k];

        int64_t bcoeff = _m_no_points - 1;
        for (; k >= 0; k--)
        {
          while (_m_binomials[bcoeff][k + 1] > idx)
          {
            bcoeff--;
          }

          indices.push_back(bcoeff);

          //std::cout << "bcoeff = " << bcoeff << std::endl;

          idx -= _m_binomials[bcoeff][k + 1];

          bcoeff--;
        }
        std::vector<int64_t>::iterator it = indices.begin(), it2;

        while (it != indices.end())
        {
          int64_t new_idx;
          conversion_with_skip(indices.begin(), indices.end(), it, new_idx);
          boundary.push_back(new_idx);
          it++;
        }
        return;
      }


      void _get_local_coboundary(int64_t idx, std::vector< int64_t >& coboundary) const
      {

        //std::cout << "glc " << idx << std::endl;

        coboundary.clear();
        int64_t dim = _get_local_dim(idx);
        //std::cout << "dim=" << dim << std::endl;
        if (dim == _m_upper_dim)
        {
          return;
        }

        static std::vector<int64_t> point_indices;
        point_indices.clear();
        conversion(idx, std::back_inserter(point_indices)); // The point indices are sorted in drcreasing order

        //std::cout << "Found " << point_indices.size() << " points" << std::endl;

        int64_t no_boundary_points = point_indices.size();
        int64_t index_to_check = point_indices.size() - 1; // used to avoid repetitions

        for (int64_t point_idx = 0; point_idx < _m_no_points; point_idx++)
        {
          if (index_to_check >= 0 && point_idx == point_indices[index_to_check])
          {
            // This point already appears in the boundary, so skip it
            index_to_check--;
            continue;
          }

          /*
          std::cout << "Instance: ";
          for(int i=0;i<point_indices.size();i++) {
          std::cout << point_indices[i] << " ";
          }
          std::cout << "Extra: " << point_idx << " ";
          */

          int64_t coboundary_idx;
          conversion_with_extra_index(point_indices.begin(), point_indices.end(), point_idx, coboundary_idx);
          //std::cout << " --> " << coboundary_idx << std::endl;
          coboundary.push_back(coboundary_idx);
        }
        /*
        std::cout << "Result: ";
        for(int i=0;i<coboundary.size();i++) {
        std::cout << coboundary[i] << " ";
        }
        std::cout << std::endl;
        */
        return;
      }

    public:
      void _get_global_dims(const std::vector< int64_t >& queries,
                            std::vector< int64_t >& answers) const
      {
        answers.resize(queries.size());
        for (int64_t idx = 0; idx < (int64_t)queries.size(); idx++)
          answers[idx] = _get_local_dim(queries[idx]);

      }

      void _get_global_boundaries(const std::vector< int64_t >& queries,
                                  data_structures::write_once_array_of_arrays< int64_t >& answers) const
      {
        std::vector< int64_t > boundary;
        for (int64_t idx = 0; idx < (int64_t)queries.size(); idx++)
        {
          _get_local_boundary(queries[idx], boundary);
          answers.set(idx, boundary.begin(), boundary.end());
        }
      }

      void _get_global_coboundaries(const std::vector< int64_t >& queries,
                                    data_structures::write_once_array_of_arrays< int64_t >& answers) const
      {
        std::vector< int64_t > coboundary;
        for (int64_t idx = 0; idx < (int64_t)queries.size(); idx++)
        {
          _get_local_coboundary(queries[idx], coboundary);
          answers.set(idx, coboundary.begin(), coboundary.end());
        }
      }


    protected:

      void _precompute_breakpoints()
      {
        _m_breakpoints.push_back(0);
        for (int64_t i = 1; i <= _m_upper_dim + 1; i++)
        {
          _m_breakpoints.push_back(_m_breakpoints[i - 1] + _m_binomials[_m_no_points][i]);
        }
        /*
        std::cout << "---" << std::endl;
        for(int i=0;i<_m_breakpoints.size();i++) {
        std::cout << _m_breakpoints[i] << std::endl;
        }
        std::cout << "---" << std::endl;
        */

      }

      void _precompute_binomials()
      {
        for (int64_t i = 0; i <= _m_no_points; i++)
        {
          std::vector<int64_t> binom_row;
          for (int64_t j = 0; j <= _m_upper_dim + 1; j++)
          {
            binom_row.push_back(binom(i, j));
          }
          _m_binomials.push_back(binom_row);
        }
      }

      // Very naive computation, todo: boost (or something else)
      int64_t binom(int64_t n, int64_t k) const
      {
        if (n < k)
        {
          return 0;
        }
        int64_t result = 1;
        for (int64_t i = n - k + 1; i <= n; i++)
        {
          result *= i;
        }
        for (int64_t i = k; i >= 1; i--)
        {
          result /= i;
        }
        return result;
      }

      // index -> set of points (in decreasing index order)
      template<typename OutputIterator>
      OutputIterator conversion(int64_t idx, OutputIterator out) const
      {
        int64_t k = _get_local_dim(idx);
        idx -= _m_breakpoints[k];
        int64_t bcoeff = _m_no_points - 1;
        for (; k >= 0; k--)
        {
          while (_m_binomials[bcoeff][k + 1] > idx)
          {
            bcoeff--;
          }
          *out++ = bcoeff;
          idx -= _m_binomials[bcoeff][k + 1];
          bcoeff--;
        }
        return out;
      }

      template<typename InputIterator>
      double diameter(InputIterator begin, InputIterator end) const
      {
        if (begin == end)
        {
          return 0.;
        }
        double max = 0.;
        InputIterator curr = begin;
        do
        {
          for (InputIterator run = curr + 1; run != end; run++)
          {
            double cdist = _m_distance_matrix[*curr][*run];
            if (cdist > max)
            {
              max = cdist;
            }
          }
          curr++;
        } while (curr != end);
        return max / 2;
      }

      // Assumes that [begin,end) is sorted in decreasing order!
      template<typename InputIterator>
      int64_t conversion_with_skip(InputIterator begin, InputIterator end, InputIterator skip, int64_t& ind) const
      {
        int64_t dist = std::distance(begin, end);
        if (skip != end)
        {
          dist--;
        }
        assert(dist <= (int)_m_upper_dim);
        InputIterator it = begin;
        ind = 0;
        for (int64_t k = dist; k >= 1; k--)
        {
          if (it == skip)
          {
            it++;
          }
          ind += _m_binomials[*it++][k];
        }
        ind += _m_breakpoints[dist - 1];
        return ind;
      }

      // Assumes that [begin,end) is sorted in decreasing order!
      template<typename InputIterator>
      int64_t conversion_with_extra_index(InputIterator begin, InputIterator end, int64_t extra, int64_t& ind) const
      {
        int64_t dist = std::distance(begin, end);
        assert(dist <= (int)_m_upper_dim);
        InputIterator it = begin;
        ind = 0;
        bool handled_extra = false;
        for (int64_t k = dist + 1; k >= 1; k--)
        {
          if (!handled_extra && (it == end || extra > *it))
          {
            ind += _m_binomials[extra][k];
            handled_extra = true;
          }
          else
          {
            ind += _m_binomials[*it++][k];
          }
        }
        ind += _m_breakpoints[dist];
        return ind;
      }

    };

  }
}

