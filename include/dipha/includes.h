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

// STL includes
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <cstdint>
#include <climits>
#include <cassert>
#include <random>

// MPI include
#include <mpi.h>

// external include for memory usage statistics
#include <mem_usage/mem_usage.h>

// external include for distributed sorting
#include <psort-1.0/src/psort.h>

// DIPHA includes
#include <dipha/globals.h>
#include <dipha/mpi_utils.h>
#include <dipha/element_distribution.h>
#include <dipha/macros.h>
#include <dipha/file_types.h>

#include <dipha/data_structures/heap_column.h>
#include <dipha/data_structures/distributed_vector.h>
#include <dipha/data_structures/flat_column_stack.h>
#include <dipha/data_structures/write_once_array_of_arrays.h>
#include <dipha/data_structures/write_once_column_array.h>

#include <dipha/inputs/abstract_weighted_cell_complex.h>
#include <dipha/inputs/weighted_cubical_complex.h>
#include <dipha/inputs/weighted_explicit_complex.h>
#include <dipha/inputs/full_rips_complex.h>
#include <dipha/inputs/sparse_rips_complex.h>

#include <dipha/algorithms/get_filtration_to_cell_map.h>
#include <dipha/algorithms/get_cell_to_filtration_map.h>
#include <dipha/algorithms/generate_unreduced_columns.h>
#include <dipha/algorithms/reduction_kernel.h>
#include <dipha/algorithms/dualize_explicit_complex.h>
#include <dipha/algorithms/compute_reduced_columns.h>

#include <dipha/outputs/save_persistence_diagram.h>
