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
  namespace globals
  {
    // MPI has some 32bit limitations, so we call it blockwise
    int MPI_BLOCK_SIZE = 1 << 23;

    /** To reduce peak memory we gather some data blockwise:
        - shared memory: use a small value to maximize cache efficiency
        - distributed memory: use a large value to hide latency of interconnect */
    int64_t DIPHA_BLOCK_SIZE = 1LL << 15;

    // accumulated total number of bytes received using the mpi_utils functions (i.e. does not include the distributed sorting)
    int64_t bytes_received = 0;

    // accumulated running time for reduction_kernel
    double reduction_kernel_running_time = 0;

    // enable benchmarking - will decrease performance a bit due to more mpi_barriers
    bool benchmark = false;
  }
}